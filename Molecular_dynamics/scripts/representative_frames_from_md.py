import os
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import MDAnalysis as mda
import matplotlib.pyplot as plt
from scipy.stats import shapiro
from sklearn.decomposition import PCA
from MDAnalysis.analysis.rms import RMSD


warnings.filterwarnings("ignore", category=DeprecationWarning, module="MDAnalysis.coordinates.DCD")
warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')
warnings.filterwarnings('ignore', category=UserWarning, module='scipy')


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
        print(f"RMSD data saved to {csv_path}")
        print("\n")
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


    def plot_rmsd_distribution(self):
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



