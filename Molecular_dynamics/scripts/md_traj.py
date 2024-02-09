import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
from sklearn.decomposition import PCA
import numpy as np
import prody as pdy


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

    def analyze_rmsd_distribution(self):
        from scipy.stats import shapiro
        shapiro_test = shapiro(self.rmsd_results['RMSD'].values)
        alpha = 0.05
        if shapiro_test[1] < alpha:
            print("Rejeitamos a hipótese nula: os dados não seguem uma distribuição normal.")
        else:
            print("Não rejeitamos a hipótese nula: os dados podem seguir uma distribuição normal.")
        estatisticas_descritivas = self.rmsd_results['RMSD'].describe()
        print(estatisticas_descritivas)
        print("\n")
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
        
        # Adicionando linhas verticais tracejadas e texto para os frames selecionados
        for _, row in self.selected_frames.iterrows():
            frame = row['frame']
            descricao = row['descricao'].replace('%', '')  # Remove o símbolo '%' da descrição
            plt.axvline(x=frame, color='red', linestyle='--', linewidth=2)
            # Adicionando descrição acima das linhas verticais
            plt.text(frame, plt.ylim()[1]*0.95, f'{descricao}', rotation=45, color='red', ha='right', va='top')

        plt.title('RMSD ao Longo da Trajetória com Pontos dos PDBs')
        plt.xlabel('Frame')
        plt.ylabel('RMSD (Å)')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.grid(True)
        plt.tight_layout()
        plt.savefig('rmsd_with_points.png')
        plt.show()


    def plot_rmsd_distribution(self):
        plt.figure(figsize=(10, 6))
        sns.histplot(self.rmsd_results['RMSD'], kde=True, bins=30, color='blue', label='Distribuição de RMSD')
        
        for _, row in self.selected_frames.iterrows():
            rmsd_value = row['rmsd']
            plt.axvline(x=rmsd_value, color='red', linestyle='--', linewidth=2)
            plt.text(rmsd_value, plt.ylim()[1]*0.95, f'{row["descricao"]}', rotation=45, color='red', ha='right')

        plt.title('Distribuição dos Valores de RMSD com PDBs Representativos')
        plt.xlabel('RMSD (Å)')
        plt.ylabel('Frequência')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.grid(True)
        plt.tight_layout()
        plt.savefig('rmsd_distribution_with_points.png')
        plt.show()

    def calculate_pca(self):
        # Selecionando os átomos para PCA (usando todos os átomos ou apenas o backbone, por exemplo)
        atoms_to_analyze = self.u.select_atoms("backbone")
        
        # Inicializando a matriz para armazenar as coordenadas
        n_frames = len(self.u.trajectory)
        n_atoms = len(atoms_to_analyze)
        coordinates = np.zeros((n_frames, n_atoms * 3))
        
        # Preenchendo a matriz com as coordenadas dos átomos em cada frame
        for i, ts in enumerate(self.u.trajectory):
            coordinates[i, :] = atoms_to_analyze.positions.flatten()
        
        # Realizando PCA com 3 componentes
        self.pca = PCA(n_components=3)
        self.pca_result = self.pca.fit_transform(coordinates)

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
        plt.savefig('pca_projections.png')
        plt.show()
        



# Demonstrando o uso da classe MDTraj
def main():
    # Inicializando a classe com o caminho da topologia e da trajetória
    md_traj = MDTraj("/media/leon/FEDF-FDB3/md_thil_10replicates_100ns/1_replica/water_remov/1/5cc8_wr_1.prmtop",
                     "/media/leon/FEDF-FDB3/md_thil_10replicates_100ns/1_replica/water_remov/1/5cc8_wr_1.dcd")
    
    # Calculando RMSD
    md_traj.calculate_rmsd()
    
    # Analisando a distribuição do RMSD
    estatisticas_descritivas = md_traj.analyze_rmsd_distribution()
    
    # Selecionando frames representativos
    md_traj.select_representative_frames(estatisticas_descritivas)
    
    # Salvando PDBs representativos
    md_traj.save_representative_pdbs()
    
    # Plotando RMSD ao longo da trajetória
    md_traj.plot_rmsd()
    
    # Plotando a distribuição dos valores de RMSD
    md_traj.plot_rmsd_distribution()

    # Calculando PCA
    md_traj.calculate_pca()
    
    # Plotando resultados da PCA
    md_traj.plot_pca_projections()

if __name__ == "__main__":
    main()
