# ressaltar que o script esta otimizado para trabalhar em paralelo para conseguir processar grandes volumes de dados
import os
import warnings
import numpy as np
import pandas as pd # realizar estas operacoes com numpy
import seaborn as sns # realizar estas operacoes com matplot
import MDAnalysis as mda
import matplotlib.pyplot as plt
from sklearn.decomposition import IncrementalPCA # realizar com numpy
from MDAnalysis.analysis.rms import RMSD
from tqdm import tqdm


warnings.filterwarnings("ignore")

class MDTraj:
    def __init__(self, topology_path, trajectory_path):
        self.u = mda.Universe(topology_path, trajectory_path)
        self.rmsd_results = None
        self.selected_frames = None



    def calculate_rmsd(self):
        backbone_atoms = self.u.select_atoms("backbone and (name CA or name C or name O or name N)")
        rmsd_analysis = RMSD(backbone_atoms, reference=backbone_atoms, ref_frame=0)
        rmsd_analysis.run()
        self.rmsd_results = pd.DataFrame(rmsd_analysis.results.rmsd, columns=["Frame", "Time (ps)", "RMSD"])
        csv_path = "./1_rmsd.csv"
        self.rmsd_results.to_csv(csv_path, index=False)
        print(f"RMSD data saved to {csv_path}\n")
        estatisticas_descritivas = self.rmsd_results['RMSD'].describe()
        print(f"Estatísticas Descritivas do RMSD:\n{estatisticas_descritivas}")
        return estatisticas_descritivas

    def select_representative_frames(self, estatisticas_descritivas):
        valores_de_interesse = {
            'min': estatisticas_descritivas['min'],
            '25': estatisticas_descritivas['25%'],
            '50': estatisticas_descritivas['50%'],
            '75': estatisticas_descritivas['75%'],
            'max': estatisticas_descritivas['max']
        }
        frames_selecionados = []
        for descricao, valor in valores_de_interesse.items():
            idx = (self.rmsd_results['RMSD'] - valor).abs().idxmin()
            frame_selecionado = self.rmsd_results.iloc[idx]
            frames_selecionados.append({
                'frame': frame_selecionado['Frame'],
                'rmsd': frame_selecionado['RMSD'],
                'descricao': descricao
            })
        self.selected_frames = pd.DataFrame(frames_selecionados)

    def save_representative_pdbs(self):
        output_directory = "pdb_md_representative"
        os.makedirs(output_directory, exist_ok=True)
        for index, row in self.selected_frames.iterrows():
            descricao = row['descricao']
            frame_index = int(row['frame'])
            self.u.trajectory[frame_index]
            pdb_filename = os.path.join(output_directory, f"{descricao}_frame_{frame_index}.pdb")
            self.u.atoms.write(pdb_filename)
        print("Arquivos PDB dos frames selecionados foram salvos com as descrições ajustadas.")

    def plot_rmsd(self):
        plt.figure(figsize=(10, 6))
        plt.plot(self.rmsd_results['Frame'], self.rmsd_results['RMSD'], label='RMSD', color='blue')
        
        # Adicionando círculos para os frames selecionados
        for _, row in self.selected_frames.iterrows():
            frame = row['frame']
            rmsd_value = self.rmsd_results['RMSD'][self.rmsd_results['Frame'] == frame].iloc[0]
            plt.scatter(frame, rmsd_value, color='red', s=50, zorder=5)  # s é o tamanho do marcador

        plt.title('RMSD ao Longo da Trajetória com Pontos dos PDBs')
        plt.xlabel('Frame')
        plt.ylabel('RMSD (Å)')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.grid(True)
        plt.tight_layout()
        plt.savefig('./rmsd_with_points.png')
        plt.show()


    def plot_rmsd_histogram(self):
        plt.figure(figsize=(10, 6))
        sns.histplot(self.rmsd_results['RMSD'], kde=True, bins=30, color='blue', label='Distribuição de RMSD')
        
        for _, row in self.selected_frames.iterrows():
            rmsd_value = row['rmsd']
            plt.scatter(rmsd_value, 0, color='red', s=50, zorder=5)  # s é o tamanho do marcador, zorder garante que o marcador fique visível acima do histograma

        plt.title('Distribuição dos Valores de RMSD com PDBs Representativos')
        plt.xlabel('RMSD (Å)')
        plt.ylabel('Frequência')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.grid(True)
        plt.tight_layout()
        plt.savefig('./rmsd_distribution_with_points.png')
        plt.show()


    def calculate_pca_incremental(self, n_components=3, chunk_size=1024):
        atoms_to_analyze = self.u.select_atoms("backbone")
        pca = IncrementalPCA(n_components=n_components)
        n_frames = len(self.u.trajectory)
        n_atoms = len(atoms_to_analyze)

        # Adicionando tqdm para barra de progresso
        for start_frame in tqdm(range(0, n_frames, chunk_size), desc="Calculating PCA"):
            end_frame = min(start_frame + chunk_size, n_frames)
            coordinates_chunk = np.zeros((end_frame - start_frame, n_atoms * 3))
            for i, ts in enumerate(self.u.trajectory[start_frame:end_frame]):
                coordinates_chunk[i, :] = atoms_to_analyze.positions.flatten()
            pca.partial_fit(coordinates_chunk)

        self.u.trajectory[0]  # Reset para início
        transformed_data = np.zeros((n_frames, n_components))
	
        # Adicionando tqdm ao loop para a transformação dos dados pelo modelo PCA ajustado
        for start_frame in tqdm(range(0, n_frames, chunk_size), desc="Transforming data with PCA"):
            end_frame = min(start_frame + chunk_size, n_frames)
            coordinates_chunk = np.zeros((end_frame - start_frame, n_atoms * 3))
            for i, ts in enumerate(self.u.trajectory[start_frame:end_frame]):
                coordinates_chunk[i, :] = atoms_to_analyze.positions.flatten()
            transformed_data[start_frame:end_frame, :] = pca.transform(coordinates_chunk)

        self.pca_result = transformed_data
        return self.pca_result
        

    
    def extract_pca_collective_variables(self):
        """
        Salva as projeções PCA em arquivos TSV.
        """
        frames = np.arange(len(self.u.trajectory)).reshape(-1, 1)
        for i in range(3):  # Para as três primeiras componentes principais
            data = np.hstack((frames, self.pca_result[:, i:i+1]))
            filename = f"cv_pca{i+1}.tsv"
            pd.DataFrame(data).to_csv(filename, sep='\t', index=False, header=False)
    
    def plot_pca_projections(self):
        # Plotando a projeção dos frames nos três primeiros componentes principais
        plt.figure(figsize=(10, 6))
        plt.plot(self.pca_result[:, 0], label='PC1', color='red')
        plt.plot(self.pca_result[:, 1], label='PC2', color='green')
        plt.plot(self.pca_result[:, 2], label='PC3', color='blue')
        plt.title('Projeção dos Frames nas Três Principais Componentes da PCA')
        plt.xlabel('Frame')
        plt.ylabel('Projeção nos Componentes Principais')
        plt.legend()
        plt.grid(True)
        plt.savefig('./pca_projections.png')
        plt.show()


        
    def main(self):
        print("Calculating RMSD...")
        rmsd_statistics = self.calculate_rmsd()
        self.select_representative_frames(rmsd_statistics)
        self.save_representative_pdbs()
        self.plot_rmsd()
        self.plot_rmsd_histogram()
        print("Calculating PCA...")
        self.calculate_pca_incremental()
        self.extract_pca_collective_variables()
        self.plot_pca_projections()
        print("Finish!")
        

def main():
    md_traj = MDTraj("/media/leon/FEDF-FDB3/md_thil_10replicates_100ns/1_replica/water_remov/traj_concatenate_aligned/5cc8_wr_1.prmtop",
                "/media/leon/FEDF-FDB3/md_thil_10replicates_100ns/1_replica/water_remov/traj_concatenate_aligned/all_traj_aligned.dcd")
    
    md_traj.main()


if __name__ == "__main__":
    main()
