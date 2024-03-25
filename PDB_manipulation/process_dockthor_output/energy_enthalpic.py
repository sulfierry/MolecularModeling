import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from mpl_toolkits.mplot3d import Axes3D  # Importa o toolkit 3D

class EnthalpicInterpolation:
    def __init__(self, cv1_path, cv2_path, energy_path):
        self.cv1_data = self.load_data(cv1_path)
        self.cv2_data = self.load_data(cv2_path)
        self.energy_values = self.load_data(energy_path, energy=True)

    def load_data(self, file_path, energy=False):
        """Carrega os dados do arquivo especificado."""
        data = np.loadtxt(file_path, usecols=[1])
        return data

    def perform_kde_interpolation(self):
        """Realiza a interpolação KDE nos dados das variáveis coletivas, ponderada pelos valores de energia."""
        weights = np.exp(-self.energy_values)
        self.kde_result = gaussian_kde(dataset=np.vstack([self.cv1_data, self.cv2_data]), weights=weights)

    def calculate_energy_surface(self, xlim, ylim):
        """Calcula a superfície de energia usando os limites especificados."""
        x = np.linspace(xlim[0], xlim[1], 100)
        y = np.linspace(ylim[0], ylim[1], 100)
        xx, yy = np.meshgrid(x, y)
        zz = np.reshape(self.kde_result(np.vstack([xx.ravel(), yy.ravel()])), xx.shape)
        return xx, yy, zz

    def plot_energy_surface_3D(self, xx, yy, zz):
        """Plota a superfície de energia 3D."""
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        # Usar a superfície diretamente não mostra bem os mínimos, então vamos usar -log(zz)
        zz_log = -np.log(zz)
        surf = ax.plot_surface(xx, yy, zz_log, cmap='viridis', linewidth=0, antialiased=False)
        fig.colorbar(surf, shrink=0.5, aspect=5, label='-log(Free Energy Density)')
        ax.set_xlabel('CV1')
        ax.set_ylabel('CV2')
        ax.set_zlabel('-log(Free Energy Density)')
        plt.title('3D Free Energy Surface with KDE Interpolation')
        plt.show()

def main():
    cv1_path = './pca_cv1.txt'
    cv2_path = './pca_cv2.txt'
    energy_path = './processed_energy_values.txt'

    interp = EnthalpicInterpolation(cv1_path, cv2_path, energy_path)
    interp.perform_kde_interpolation()

    xlim = (-8, 8)  # Limites para CV1
    ylim = (-5, 7)  # Limites para CV2

    xx, yy, zz = interp.calculate_energy_surface(xlim, ylim)

    interp.plot_energy_surface_3D(xx, yy, zz)

if __name__ == "__main__":
    main()
