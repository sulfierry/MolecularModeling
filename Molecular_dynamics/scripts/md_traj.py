import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings

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

