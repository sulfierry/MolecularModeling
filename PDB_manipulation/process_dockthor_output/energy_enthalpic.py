import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from mpl_toolkits.mplot3d import Axes3D

class EnthalpicInterpolation:
    def __init__(self, cv1_path, cv2_path, energy_path):
        self.cv1_data = np.loadtxt(cv1_path, usecols=[1])
        self.cv2_data = np.loadtxt(cv2_path, usecols=[1])
        self.energy_values = np.loadtxt(energy_path, usecols=[1])

    def perform_kde_interpolation(self, bandwidth=None):
        """Realiza a interpolação KDE nos dados das variáveis coletivas, ponderada pelos valores de energia."""
        weights = np.exp(-self.energy_values)
        self.kde_result = gaussian_kde(dataset=np.vstack([self.cv1_data, self.cv2_data]), weights=weights, bw_method=bandwidth)

    def calculate_energy_surface(self, xlim, ylim):
        """Calcula a superfície de energia usando os limites especificados e inverte os valores."""
        x = np.linspace(xlim[0], xlim[1], 100)
        y = np.linspace(ylim[0], ylim[1], 100)
        xx, yy = np.meshgrid(x, y)
        zz = -np.reshape(self.kde_result(np.vstack([xx.ravel(), yy.ravel()])), xx.shape)
        return xx, yy, zz

    def plot_energy_surface_2D(self, xx, yy, zz, n_min_points=0):
        """Plota a superfície de energia 2D."""
        plt.figure(figsize=(10, 8))
        cp = plt.contourf(xx, yy, zz, levels=50, cmap='viridis', extend='both')
        plt.colorbar(cp, label='Enthalpic Density Value')
        plt.xlabel('CV1')
        plt.ylabel('CV2')
        plt.title('2D Enthalpic Energy Surface with KDE Interpolation')

        # Destaca os n_min_points com menores valores de entalpia, se necessário
        if n_min_points > 0:
            self.highlight_points(n_min_points, '2D', plt.gca())

        plt.show()


    def plot_energy_surface_3D(self, xx, yy, zz, n_min_points=None):
        """Plota a superfície de energia 3D."""
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Plota a superfície de densidade entálpica
        surf = ax.plot_surface(xx, yy, zz, cmap='viridis', linewidth=0, antialiased=False, alpha=0.8)
        fig.colorbar(surf, shrink=0.5, aspect=5, label='Enthalpic Value')
        ax.set_xlabel('CV1')
        ax.set_ylabel('CV2')
        ax.set_zlabel('Enthalpic Density Value')
        plt.title('3D Enthalpic Energy Surface with KDE Interpolation')

        # Destaca os n_min_points com menores valores de entalpia, se necessário
        if n_min_points > 0:
            # Ajuste aqui para garantir que os pontos sejam plotados de forma que não interfiram com a superfície
            # Por exemplo, você pode ajustar a cor, tamanho ou marcador dos pontos
            self.highlight_points(n_min_points, '3D', ax)

        # Ajusta a visão para garantir que a superfície seja visível
        ax.view_init(elev=20, azim=45)

        plt.show()




    def highlight_points(self, n_min_points, plot_type, ax=None):
        """Destaca os n_min_points de menor valor entálpico."""
        sorted_indices = np.argsort(self.energy_values)[:n_min_points]
        colors = ['red', 'green', 'blue', 'purple', 'orange']
        markers = ['o', '^', 's', 'P', '*']

        for i, idx in enumerate(sorted_indices):
            if plot_type == '2D':
                plt.scatter(self.cv1_data[idx], self.cv2_data[idx], color=colors[i % len(colors)],
                            marker=markers[i % len(markers)], s=100, edgecolors='black', label=f'Point {i+1}')
            elif plot_type == '3D' and ax is not None:
                ax.scatter(self.cv1_data[idx], self.cv2_data[idx], self.energy_values[idx], color=colors[i % len(colors)],
                           marker=markers[i % len(markers)], s=50, label=f'Point {i+1}')

        if plot_type == '2D':
            plt.legend()

    def calculate_density_at_point(self, cv1, cv2):
        """Calcula a densidade no ponto (cv1, cv2) usando a interpolação KDE."""
        # Supondo que self.kde_result já foi calculado
        density = self.kde_result.evaluate([cv1, cv2])
        return density


    def save_to_csv(self, output_path):
        """Salva os dados em um arquivo .csv."""
        # Cria um DataFrame com os dados
        data = pd.DataFrame({
            'Frame': np.arange(len(self.cv1_data)),
            'cv1': self.cv1_data,
            'cv2': self.cv2_data,
            'entalpia': self.energy_values,
            'densidade': [self.calculate_density_at_point(cv1, cv2) for cv1, cv2 in zip(self.cv1_data, self.cv2_data)]
        })

        # Salva o DataFrame em um arquivo .csv
        data.to_csv(output_path, index=False)




def main(n_min_points_desired):
    cv1_path = './pca_cv1.txt'
    cv2_path = './pca_cv2.txt'
    energy_path = './processed_energy_values.txt'

    interp = EnthalpicInterpolation(cv1_path, cv2_path, energy_path)
    bandwidth = 0.25  # Exemplo de valor; ajuste conforme necessário
    interp.perform_kde_interpolation(bandwidth=bandwidth)

    xlim = (-3.8, 5.5)
    ylim = (-2.2, 3)

    xx, yy, zz = interp.calculate_energy_surface(xlim, ylim)

    interp.plot_energy_surface_2D(xx, yy, zz, n_min_points_desired)
    interp.plot_energy_surface_3D(xx, yy, zz, n_min_points_desired)

    # Salva os dados em um arquivo .csv
    output_path = './output.csv'
    interp.save_to_csv(output_path)

if __name__ == "__main__":
    main(0)  # Define o número desejado de pontos mínimos a serem destacados
