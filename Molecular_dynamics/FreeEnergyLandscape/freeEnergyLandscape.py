import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from matplotlib.colors import LinearSegmentedColormap



class FreeEnergyLandscape:
    
    def __init__(self, cv1_path, cv2_path, temperature=300, boltzmann_constant=8.314e-3):
        self.cv1_path = cv1_path
        self.cv2_path = cv2_path
        self.temperature = temperature
        self.kB = boltzmann_constant
        self.colors = [
            (0, "darkblue"),    # 0 a 3
            (3/25, "blue"),     # 3 a 6
            (6/25, "lightblue"),# 6 a 9
            (9/25, "#ADD8E6"),  # 9 a 12 azul claríssimo
            (12/25, "#FFA07A"), # 12 a 15 vermelho claro (quase salmão)
            (15/25, "#FF4500"), # 15 a 18 mais escuro (quase laranja)
            (18/25, "#FF6347"), # 18 a 21 laranja/vermelho
            (21/25, "darkred"), # 21 a 24 vermelho escuro
            (1, "darkred")      # Garante que o máximo seja vermelho escuro
        ]
        self.custom_cmap = LinearSegmentedColormap.from_list("custom_energy", self.colors)
        self.proj1_data_original = None
        self.proj2_data_original = None


    def load_data(self):
        self.proj1_data_original = np.loadtxt(self.cv1_path, usecols=[1])
        self.proj2_data_original = np.loadtxt(self.cv2_path, usecols=[1])
