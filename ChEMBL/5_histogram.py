import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from scipy.spatial.distance import cosine
import seaborn as sns

class Histogram:
	"""
	O presente algoritmo compara todos os fingerprints moleculares uns contra os outros, 
	o que é conhecido como uma abordagem "todos contra todos". Para cada fingerprint no 
	conjunto de dados,ele calcula a similaridade ou a distância com todos os outros 
	fingerprints. Este processo é repetido para cada fingerprint, resultando em um 
	conjunto abrangente de medidas de similaridade e distância entre cada par de moléculas 
	representadas pelos fingerprints.

	Essa abordagem pode ser computacionalmente intensiva, especialmente para conjuntos de dados 
	grandes, porque o número de comparações cresce quadraticamente com o número de fingerprints. 
	Por exemplo, se há N fingerprints, o número de comparações será N×(N−1)/2, o que pode ser um 
	número muito grande para grandes conjuntos de dados. É por isso que o algoritmo processa os 
	dados em lotes, para gerenciar o uso de memória e recursos de processamento.
	"""

    def __init__(self, file_path, batch_size=1024):
        self.file_path = file_path
        self.batch_size = batch_size

    @staticmethod
    def smiles_to_fingerprint(smiles, radius=2):
        mol = Chem.MolFromSmiles(smiles)
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius) if mol else None

    @staticmethod
    def process_batch(fingerprints, calculate_function):
        results = []
        for fp in fingerprints:
            results.extend(calculate_function(fp, fingerprints))
        return results

    @staticmethod
    def calculate_similarity(fingerprint, fingerprints, similarity_metric):
        similarities = []
        for fp in fingerprints:
            if similarity_metric == 'cosine':
                arr1, arr2 = np.zeros((1,)), np.zeros((1,))
                DataStructs.ConvertToNumpyArray(fingerprint, arr1)
                DataStructs.ConvertToNumpyArray(fp, arr2)
                similarity = 1 - cosine(arr1, arr2)
            elif similarity_metric == 'tanimoto':
                similarity = DataStructs.TanimotoSimilarity(fingerprint, fp)
            elif similarity_metric == 'dice':
                similarity = DataStructs.DiceSimilarity(fingerprint, fp)
            similarities.append(similarity)
        return similarities

    @staticmethod
    def calculate_distance(fingerprint, fingerprints, distance_metric):
        distances = []
        for fp in fingerprints:
            if distance_metric == 'hamming':
                dist = sum(1 for a, b in zip(fingerprint, fp) if a != b)
            elif distance_metric == 'manhattan':
                arr1, arr2 = np.zeros((1,)), np.zeros((1,))
                DataStructs.ConvertToNumpyArray(fingerprint, arr1)
                DataStructs.ConvertToNumpyArray(fp, arr2)
                dist = np.sum(np.abs(arr1 - arr2))
            distances.append(dist)
        return distances

    def process_in_batches(self, data, process_function, metric):
        results = []
        for i in range(0, len(data), self.batch_size):
            batch = data[i:i+self.batch_size]
            results.extend(self.process_batch(batch, lambda fp, fps: process_function(fp, fps, metric)))
        return results

    def create_histograms(self):
        data = pd.read_csv(self.file_path, sep='\t')
        data['fingerprints'] = data['canonical_smiles'].apply(self.smiles_to_fingerprint)
        valid_fps = data['fingerprints'].dropna().tolist()

        similarity_metrics = ['tanimoto', 'dice', 'cosine']
        distance_metrics = ['hamming', 'manhattan']
        all_similarities = {metric: self.process_in_batches(valid_fps, self.calculate_similarity, metric) for metric in similarity_metrics}
        all_distances = {metric: self.process_in_batches(valid_fps, self.calculate_distance, metric) for metric in distance_metrics}

        fig, axs = plt.subplots(3, 2, figsize=(13, 13))
        for ax, (metric, values) in zip(axs[:, 0], all_similarities.items()):
            ax.hist(values, bins=50, color='skyblue', edgecolor='black')
            ax.set_title(f'Similarity - {metric.capitalize()}', fontsize=10)
            ax.set_ylabel('Frequency')

        for ax, (metric, values) in zip(axs[:, 1], all_distances.items()):
            ax.hist(values, bins=50, color='salmon', edgecolor='black')
            ax.set_title(f'Distance - {metric.capitalize()}', fontsize=10)
            ax.set_ylabel('Frequency')

        axs[2, 1].hist(data['pchembl_value'].dropna(), bins=20, color='green', edgecolor='black')
        axs[2, 1].set_title('ChEMBL value', fontsize=10)
        axs[2, 1].set_ylabel('Frequency')

        plt.tight_layout()
        plt.savefig('histogram_similarity_distance.png')
        plt.show()

if __name__ == "__main__":
    histogram = Histogram('./chembl_cluster_hits.tsv')
    histogram.create_histograms()

