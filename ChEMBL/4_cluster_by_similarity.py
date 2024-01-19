"""
O presente algoritimo agurpa moléculas com base em suas similaridades estruturais. 
Utilizando a representação de SMILES (Simplified Molecular Input Line Entry System) 
para moléculas, o algoritmo converte estas estruturas em fingerprints moleculares e 
realiza a clusterização baseada em um limiar de similaridade predefinido.
"""
import os
import sys
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
import matplotlib.pyplot as plt
from multiprocessing import Pool
from sklearn.manifold import TSNE
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, DataStructs

class MoleculeClusterer:
    
    def __init__(self, smiles_file_path):
        self.smiles_file_path = smiles_file_path
        self.data = None
        self.fingerprints = []

    def load_data(self, smile_column):
        # Correção aqui: usar 'smiles_file_path' em vez de 'data_path'
        self.data = pd.read_csv(self.smiles_file_path, sep='\t')
        self.data = self.data.dropna(subset=[smile_column])
    
    def cluster_by_similarity(self, threshold=0.8):
        num_fps = len(self.fingerprints)
        clusters = []
        visited = set()

        for i in range(num_fps):
            if i in visited:
                continue

            cluster = [i]
            for j in range(i + 1, num_fps):
                similarity = DataStructs.TanimotoSimilarity(self.fingerprints[i], self.fingerprints[j])
                if similarity >= threshold:
                    cluster.append(j)
                    visited.add(j)

            clusters.append(cluster)

        return clusters

    @staticmethod
    def smiles_to_fingerprint(smiles):
        mol = Chem.MolFromSmiles(smiles)
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2) if mol else None

    @staticmethod
    def compute_fingerprint(smiles):
        return MoleculeClusterer.smiles_to_fingerprint(smiles)

    def parallel_generate_fingerprints(self, smile_column):
        num_cpus = os.cpu_count()
        with Pool(num_cpus) as pool:
            self.fingerprints = list(tqdm(pool.imap(self.compute_fingerprint, self.data[smile_column]), total=len(self.data)))
        self.fingerprints = [fp for fp in self.fingerprints if fp is not None]

    def save_clustered_data(self, output_file_path):
        self.data.to_csv(output_file_path, sep='\t', index=False)

    def save_state(self, file_path):
        state = {
            'data': self.data,
            'fingerprints': self.fingerprints
        }
        with open(file_path, 'wb') as file:
            pickle.dump(state, file)

    def load_state(self, file_path):
        with open(file_path, 'rb') as file:
            state = pickle.load(file)
            self.data = state['data']
            self.fingerprints = state['fingerprints']


    def calculate_tsne(self):
        # Convertendo a lista de fingerprints para um array NumPy com uma barra de progresso do tqdm
        fingerprint_array = np.array([fp for fp in tqdm(self.fingerprints, desc='Processing Fingerprints') if fp is not None])
	tsne = TSNE(n_components=2, random_state=42, learning_rate=200.0, init='pca')

        tsne_results = tsne.fit_transform(fingerprint_array)
        return tsne_results


    def plot_tsne(self, tsne_results, threshold):
        plt.figure(figsize=(12, 8))

        # Filtrar os dados para incluir apenas clusters com tamanho >= threshold
        filtered_data = self.data[self.data['cluster_id'].map(self.data['cluster_id'].value_counts()) >= threshold]

        # Filtrar os resultados t-SNE para incluir apenas os dados filtrados
        filtered_indices = filtered_data.index
        filtered_tsne_results = tsne_results[filtered_indices]

        # Obter os tamanhos dos clusters para os dados filtrados
        cluster_sizes = filtered_data['cluster_id'].map(filtered_data['cluster_id'].value_counts())

        # Mapear os tamanhos dos clusters para um intervalo de cores
        cmap = plt.cm.viridis
        norm = plt.Normalize(cluster_sizes.min(), cluster_sizes.max())

        scatter = plt.scatter(filtered_tsne_results[:, 0], filtered_tsne_results[:, 1], 
                            c=cluster_sizes, cmap=cmap, norm=norm, alpha=0.5)

        cbar = plt.colorbar(scatter)
        cbar.set_label('Cluster size')

        plt.title('Molecular similarity in 2D space')
        plt.xlabel('t-SNE feature 1')
        plt.ylabel('t-SNE feature 2')
        plt.savefig('tsne_similarity.png')
        plt.show()
        
        
    def plot_cluster_size_distribution(self, threshold):
        # Filtrar os dados para incluir apenas clusters com tamanho >= threshold
        cluster_sizes = self.data['cluster_id'].value_counts()
        filtered_cluster_sizes = cluster_sizes[cluster_sizes >= threshold]

        plt.figure(figsize=(10, 6))
        plt.hist(filtered_cluster_sizes, bins=range(threshold, filtered_cluster_sizes.max() + 1), 
                alpha=0.7, color='blue', edgecolor='black')

        plt.title('Smiles distribution per cluster')
        plt.xlabel('Number of smiles per cluster')
        plt.ylabel('Cluster count')
        plt.grid(axis='y', alpha=0.75)
        plt.savefig('cluster_size_distribution.png')
        plt.show()


    def save_clusters_as_tsv(self, threshold, smile_column, cluster_hits_file):
        # Criar pasta 'clusters' se não existir
        os.makedirs('./clusters', exist_ok=True)

        # Calcula o peso molecular para cada SMILES e adiciona como uma nova coluna
        self.data['molecular_weight'] = self.data[smile_column].apply(lambda x: Descriptors.MolWt(Chem.MolFromSmiles(x)) if x else None)

        # Colunas a serem incluídas nos arquivos de saída
        output_columns = ['molregno', 'target_kinase', 'canonical_smiles', 'standard_value', 'standard_type', 'pchembl_value',
                        'compound_name', 'molecular_weight', 'cluster_id']
   

        # Lista para armazenar os hits de cada cluster
        cluster_hits = []

        for cluster_id in set(self.data['cluster_id']):
            cluster_data = self.data[self.data['cluster_id'] == cluster_id]

            # Renomear as colunas conforme necessário
            #cluster_data = cluster_data.rename(columns=rename_columns)

            # Ordenar os dados do cluster por peso molecular
            cluster_data_sorted = cluster_data.sort_values('molecular_weight')

            # Salvar apenas se o cluster contiver pelo menos 'threshold' SMILES
            if len(cluster_data_sorted) >= threshold:
                # Selecionar apenas as colunas desejadas para salvar
                cluster_data_to_save = cluster_data_sorted[output_columns]
                cluster_data_to_save.to_csv(f'./clusters/cluster_{cluster_id}.tsv', sep='\t', index=False)

                # Adicionar o hit do cluster (molecule com o menor peso molecular) à lista de hits
                cluster_hits.append(cluster_data_sorted.iloc[0])

        # Criar DataFrame dos hits dos clusters e salvar
        cluster_hits_df = pd.DataFrame(cluster_hits, columns=output_columns)
        cluster_hits_df.to_csv(cluster_hits_file, sep='\t', index=False)


def run(smiles_file_path, output_file_path, state_file_path, tanimoto_threshold, cluster_size_threshold, smile_column):
    
    clusterer = MoleculeClusterer(smiles_file_path)

    try:
        clusterer.load_state(state_file_path)
        print("Estado carregado com sucesso.")
    except FileNotFoundError:
        print("Nenhum estado salvo encontrado. Iniciando processamento do zero.")
        clusterer.load_data(smile_column)
        clusterer.parallel_generate_fingerprints(smile_column)
        clusterer.save_state(state_file_path)
    except Exception as e:
        print(f"Erro ao carregar o estado: {e}")
        return

    if clusterer.fingerprints:
        tsne_results = clusterer.calculate_tsne()
        clusters = clusterer.cluster_by_similarity(threshold=tanimoto_threshold)

        cluster_ids = [None] * len(clusterer.data)
        for cluster_id, cluster in enumerate(clusters):
            for idx in cluster:
                cluster_ids[idx] = cluster_id
        clusterer.data['cluster_id'] = cluster_ids

        clusterer.save_clusters_as_tsv(cluster_size_threshold, smile_column, './chembl_cluster_hits.tsv')
        clusterer.plot_tsne(tsne_results, cluster_size_threshold)
        clusterer.plot_cluster_size_distribution(cluster_size_threshold)

        #clusterer.save_clustered_data(output_file_path)
        clusterer.save_state(state_file_path)
    else:
        print("Nenhum fingerprint válido foi encontrado.")


def main():
    smiles_file_path = '../../../1_database/nr_kinase_drug_info_all_manually_validated_IC50_Ki_kd_10uM.tsv'
    output_file_path = './clustered_smiles.tsv'
    state_file_path  = './molecule_clusterer_state.pkl'

    smile_column = 'canonical_smiles'
    tanimoto_threshold = 0.8
    cluster_size_threshold = 3 # numero minimo de moleculas por cluster
    
    run(smiles_file_path, output_file_path, state_file_path, tanimoto_threshold, cluster_size_threshold, smile_column)
    
if __name__ == "__main__":
    main()
