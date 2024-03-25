import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

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

    def plot_energy_surface(self, xx, yy, zz, left, right, top, bottom):
        """Plota a superfície de energia 2D, ajustando os limites conforme especificado."""
        plt.figure(figsize=(10, 8))
        # Normalização para escala de -6 a 6 para a régua
        zz_normalized = np.interp(zz, (zz.min(), zz.max()), (-6, 6))
        cp = plt.contourf(xx, yy, zz_normalized, levels=np.linspace(-6, 6, 50), cmap='viridis', extend='both')
        plt.colorbar(cp, label='Free Energy Density')
        plt.xlim(left=left, right=right)
        plt.ylim(bottom=bottom, top=top)
        plt.xlabel('CV1')
        plt.ylabel('CV2')
        plt.title('2D Free Energy Surface with KDE Interpolation')
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

    # Limites para ajuste de visualização
    left, right, top, bottom = -3, 5, 3, -3

    interp.plot_energy_surface(xx, yy, zz, left, right, top, bottom)

if __name__ == "__main__":
    main()
