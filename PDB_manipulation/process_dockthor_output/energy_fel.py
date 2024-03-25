import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

class EnthalpicInterpolation:
    def __init__(self, cv1_path, cv2_path, energy_path):
        self.cv1_data = self.load_data(cv1_path)
        self.cv2_data = self.load_data(cv2_path)
        self.energy_values = self.load_data(energy_path, energy=True)

    def load_data(self, file_path, energy=False):
        """
        Carrega os dados do arquivo especificado.
        Se energy=True, processa o arquivo de valores de energia.
        """
        data = np.loadtxt(file_path, usecols=[1])
        return data

    def perform_kde_interpolation(self, weights=True):
        """
        Realiza a interpolação KDE nos dados das variáveis coletivas,
        opcionalmente ponderada pelos valores de energia.
        """
        if weights:
            weights = np.exp(-self.energy_values)
        else:
            weights = None
        self.kde_result = gaussian_kde(dataset=np.vstack([self.cv1_data, self.cv2_data]), weights=weights)

    def calculate_scaled_energy_surface(self, xlim, ylim, levels):
        """
        Calcula a superfície de energia ajustada e retorna os valores para plotagem.
        """
        xx, yy = np.mgrid[xlim[0]:xlim[1]:100j, ylim[0]:ylim[1]:100j]
        zz = self.kde_result(np.vstack([xx.flatten(), yy.flatten()]))
        zz = np.reshape(zz, xx.shape)

        # -log transform and scale
        zz_log = -np.log(zz)
        zz_normalized = (zz_log - zz_log.min()) / (zz_log.max() - zz_log.min())
        zz_scaled = zz_normalized * (levels[-1] - levels[0]) + levels[0]
        return xx, yy, zz_scaled

    def plot_energy_surface(self, xx, yy, zz, levels):
        """
        Plota a superfície de energia 2D.
        """
        plt.figure(figsize=(10, 8))
        cp = plt.contourf(xx, yy, zz, levels=levels, extend='both', cmap='viridis')
        plt.colorbar(cp, label='Scaled Free Energy')
        plt.xlabel('CV1')
        plt.ylabel('CV2')
        plt.title('2D Free Energy Surface with KDE Interpolation')
        plt.show()

def main():
    # Caminhos para os arquivos
    cv1_path = './pca_cv1.txt'
    cv2_path = './pca_cv2.txt'
    energy_path = './processed_energy_values.txt'

    # Criação e utilização da classe
    interp = EnthalpicInterpolation(cv1_path, cv2_path, energy_path)
    interp.perform_kde_interpolation()

    # Define limites e níveis para a escala da energia
    xlim = (interp.cv1_data.min(), interp.cv1_data.max())
    ylim = (interp.cv2_data.min(), interp.cv2_data.max())
    levels = np.arange(-6, 8, 2)  # De -6 a 6, variação de 2

    xx, yy, zz_scaled = interp.calculate_scaled_energy_surface(xlim, ylim, levels)
    interp.plot_energy_surface(xx, yy, zz_scaled, levels)

if __name__ == "__main__":
    main()
