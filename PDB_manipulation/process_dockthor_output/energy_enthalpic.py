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
        """Performs KDE interpolation on the collective variable data, weighted by the energy values."""
        weights = np.exp(-self.energy_values)
        self.kde_result = gaussian_kde(dataset=np.vstack([self.cv1_data, self.cv2_data]), weights=weights, bw_method=bandwidth)

    def calculate_energy_surface(self, xlim, ylim):
        """Calculates the energy surface using the specified limits and inverts the values."""
        x = np.linspace(xlim[0], xlim[1], 100)
        y = np.linspace(ylim[0], ylim[1], 100)
        xx, yy = np.meshgrid(x, y)
        zz = -np.reshape(self.kde_result(np.vstack([xx.ravel(), yy.ravel()])), xx.shape)
        return xx, yy, zz

    def plot_energy_surface_2D(self, xx, yy, zz, n_min_points=0):
        """Plots the 2D energy surface."""
        plt.figure(figsize=(10, 8))
        cp = plt.contourf(xx, yy, zz, levels=50, cmap='viridis', extend='both')
        plt.colorbar(cp, label='Enthalpic Density Value')
        plt.xlabel('CV1')
        plt.ylabel('CV2')
        plt.title('2D Enthalpic Energy Surface with KDE Interpolation')
        if n_min_points > 0:
            self.highlight_points(n_min_points, '2D', plt.gca())
        plt.show()

    def plot_energy_surface_3D(self, xx, yy, zz, n_min_points=None):
        """Plots the 3D energy surface."""
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(xx, yy, zz, cmap='viridis', linewidth=0, antialiased=False, alpha=0.8)
        fig.colorbar(surf, shrink=0.5, aspect=5, label='Enthalpic Value')
        ax.set_xlabel('CV1')
        ax.set_ylabel('CV2')
        ax.set_zlabel('Enthalpic Density Value')
        plt.title('3D Enthalpic Energy Surface with KDE Interpolation')
        if n_min_points > 0:
            self.highlight_points(n_min_points, '3D', ax)
        ax.view_init(elev=20, azim=45)
        plt.show()

    def highlight_points(self, n_min_points, plot_type, ax=None):
        """Highlights the n_min_points of lowest enthalpic value."""
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
        """Calculates the density at point (cv1, cv2) using KDE interpolation."""
        density = self.kde_result.evaluate([cv1, cv2])
        return density

    def save_to_csv(self, output_path):
        """Saves the data to a .csv file, with 'densidade' and 'entalpia' columns adjusted based on KDE density."""
        # Estimating density for each point
        densities = [self.calculate_density_at_point(cv1, cv2) for cv1, cv2 in zip(self.cv1_data, self.cv2_data)]

        # Creating a DataFrame with the adjusted data
        data = pd.DataFrame({
            'Frame': np.arange(len(self.cv1_data)),
            'cv1': self.cv1_data,
            'cv2': self.cv2_data,
            'densidade': densities,
            'entalpia': self.energy_values
        })

        # Sorting the DataFrame by 'densidade' to adjust 'entalpia' accordingly
        data_sorted_by_density = data.sort_values(by='densidade', ascending=False)
        entalpia_sorted = data_sorted_by_density['entalpia'].values

        # Reassigning 'entalpia' based on sorted densities
        data['entalpia'] = entalpia_sorted

        # Saving the DataFrame to a .csv file
        data.to_csv(output_path, index=False)

def main(n_min_points_desired):
    cv1_path = './pca_cv1.txt'
    cv2_path = './pca_cv2.txt'
    energy_path = './processed_energy_values.txt'

    interp = EnthalpicInterpolation(cv1_path, cv2_path, energy_path)
    bandwidth = 0.25  # Example bandwidth; adjust as needed
    interp.perform_kde_interpolation(bandwidth=bandwidth)

    xlim = (-3.8, 5.5)
    ylim = (-2.2, 3)

    xx, yy, zz = interp.calculate_energy_surface(xlim, ylim)

    interp.plot_energy_surface_2D(xx, yy, zz, n_min_points_desired)
    interp.plot_energy_surface_3D(xx, yy, zz, n_min_points_desired)

    # Saving the data to a .csv file
    output_path = './output.csv'
    interp.save_to_csv(output_path)

if __name__ == "__main__":
    main(0)  # Desired number of minimum points to be highlighted
