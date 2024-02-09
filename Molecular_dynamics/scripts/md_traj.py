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

