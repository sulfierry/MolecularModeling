import os
import sys
import shutil
import tempfile
import numpy as np
import imageio.v2 as imageio
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
# from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap

class FreeEnergyLandscape:
    
    def __init__(self, cv1_path, cv2_path, temperature, boltzmann_constant):
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

    def boltzmann_inversion_original(self, data, title):
        hist, bin_edges = np.histogram(data, bins=100, density=True)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        hist = np.clip(hist, a_min=1e-10, a_max=None)
        G = -self.kB * self.temperature * np.log(hist)
        G = np.clip(G - np.min(G), 0, 25)
        plt.figure(figsize=(10, 6))
        plt.plot(bin_centers, G, label='Energia Livre', color='red')
        plt.title(f'Paisagem Energética Livre de {title}')
        plt.xlabel('Valor')
        plt.ylabel('Energia Livre (kJ/mol)')
        plt.ylim(0, 25)
        plt.grid(True)
        plt.legend()
        plt.show()

    def plot_energy_landscape(self):
        values_original = np.vstack([self.proj1_data_original, self.proj2_data_original])
        kernel_original = gaussian_kde(values_original)
        X_original, Y_original = np.mgrid[self.proj1_data_original.min():self.proj1_data_original.max():100j, 
                                          self.proj2_data_original.min():self.proj2_data_original.max():100j]
        positions_original = np.vstack([X_original.ravel(), Y_original.ravel()])
        Z_original = np.reshape(kernel_original(positions_original).T, X_original.shape)
        G_original = -self.kB * self.temperature * np.log(Z_original)
        G_original = np.clip(G_original - np.min(G_original), 0, 25)
        plt.figure(figsize=(8, 6))
        plt.contourf(X_original, Y_original, G_original, levels=np.linspace(0, 25, 100), cmap=self.custom_cmap)
        plt.colorbar(label='Energia Livre (kJ/mol)', ticks=range(0, 26, 3))
        plt.xlabel('CV1 (Ângulo)')
        plt.ylabel('CV2 (Distância)')
        plt.title('Paisagem Energética Gerada')
        plt.show()

    def plot_3D_energy_landscape(self):
        values_original = np.vstack([self.proj1_data_original, self.proj2_data_original])
        kernel_original = gaussian_kde(values_original)
        X_original, Y_original = np.mgrid[self.proj1_data_original.min():self.proj1_data_original.max():100j, 
                                          self.proj2_data_original.min():self.proj2_data_original.max():100j]
        positions_original = np.vstack([X_original.ravel(), Y_original.ravel()])
        Z_original = np.reshape(kernel_original(positions_original).T, X_original.shape)
        G_original = -self.kB * self.temperature * np.log(Z_original)
        G_original = np.clip(G_original - np.min(G_original), 0, 25)
        
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(X_original, Y_original, G_original, cmap=self.custom_cmap, edgecolor='none')
        ax.set_xlabel('CV1 (Ângulo)')
        ax.set_ylabel('CV2 (Distância)')
        ax.set_zlabel('Energia Livre (kJ/mol)')
        ax.set_title('Paisagem Energética em 3D')
        fig.colorbar(surf, shrink=0.5, aspect=5, label='Energia Livre (kJ/mol)')
        plt.show()

    def create_3D_gif(self, gif_filename='energy_landscape_3D.gif', n_angles=36, elevation=30):
        temp_dir = tempfile.mkdtemp()  # Cria um diretório temporário
        filenames = []

        values_original = np.vstack([self.proj1_data_original, self.proj2_data_original])
        kernel_original = gaussian_kde(values_original)
        X_original, Y_original = np.mgrid[self.proj1_data_original.min():self.proj1_data_original.max():100j, 
                                          self.proj2_data_original.min():self.proj2_data_original.max():100j]
        positions_original = np.vstack([X_original.ravel(), Y_original.ravel()])
        Z_original = np.reshape(kernel_original(positions_original).T, X_original.shape)
        G_original = -self.kB * self.temperature * np.log(Z_original)
        G_original = np.clip(G_original - np.min(G_original), 0, 25)

        for angle in range(0, 360, int(360 / n_angles)):
            fig = plt.figure(figsize=(10, 7))
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_surface(X_original, Y_original, G_original, cmap=self.custom_cmap, edgecolor='none')
            ax.view_init(elev=elevation, azim=angle)
            ax.set_xlabel('CV1 (Ângulo)')
            ax.set_ylabel('CV2 (Distância)')
            ax.set_zlabel('Energia Livre (kJ/mol)')

            # Salva a figura diretamente em um arquivo
            frame_filename = os.path.join(temp_dir, f"frame_{angle:03d}.png")
            plt.savefig(frame_filename)
            filenames.append(frame_filename)
            plt.close()

        # Cria o GIF a partir dos arquivos salvos
        with imageio.get_writer(gif_filename, mode='I', duration=0.1) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)

        # Limpa os arquivos temporários
        shutil.rmtree(temp_dir)


    def main(self):
        self.load_data()
        self.boltzmann_inversion_original(self.proj1_data_original, 'CV1 (Ângulo)')
        self.boltzmann_inversion_original(self.proj2_data_original, 'CV2 (Distância)')
        self.plot_energy_landscape()
        self.plot_3D_energy_landscape()
        self.create_3D_gif()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python freeEnergyLandscape.py path/to/cv1_data.txt path/to/cv2_data.txt")
        print("Example: python freeEnergyLandscape.py cv1_angle.txt cv2_distance.txt")
        sys.exit(1)

    cv1_path = sys.argv[1]
    cv2_path = sys.argv[2]

    try:
        t = 300         # Temperature in K
        kB = 8.314e-3   # Boltzmann constant in kJ/(mol·K)

        fel = FreeEnergyLandscape(cv1_path, cv2_path, t, kB)
        fel.main()
    except FileNotFoundError as e:
        print(f"Error: File not found. {e}")
        print("Please ensure the file paths are correct and try again.")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)
