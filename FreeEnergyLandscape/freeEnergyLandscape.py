import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from matplotlib.colors import LinearSegmentedColormap

class GaussianKDE:
    def __init__(self, dataset, bandwidth='scott'):
        self.dataset = np.atleast_2d(dataset)
        self.n, self.d = self.dataset.shape
        if bandwidth == 'scott':
            self.bandwidth = np.power(self.n, -1./(self.d+4))  # Scott's Rule
        else:
            self.bandwidth = bandwidth
        self.factor = np.sqrt(2 * np.pi) * self.bandwidth ** self.d

    def _kernel_function(self, dist_squared):
        return np.exp(-dist_squared / (2 * self.bandwidth ** 2))

    def evaluate(self, points):
        points = np.atleast_2d(points)
        n_points, d_points = points.shape
        if d_points != self.d:
            raise ValueError("Os pontos devem ter a mesma dimensão dos dados do dataset.")
        
        # Calcula as distâncias de forma mais eficiente
        dist_squared = cdist(points, self.dataset, 'sqeuclidean')

        # Calcula os valores do kernel para cada par ponto-dado
        kernel_values = self._kernel_function(dist_squared)

        # Soma as contribuições dos kernels e normaliza
        density = np.sum(kernel_values, axis=1) / (self.n * self.factor)
        return density
    
    def __call__(self, points):
        return self.evaluate(points)


class FreeEnergyLandscape:
    def __init__(self, cv1_path, cv2_path, temperature, boltzmann_constant):
        self.cv1_path = cv1_path
        self.cv2_path = cv2_path
        self.temperature = temperature
        self.kB = boltzmann_constant
        self.colors = [
            (0, "darkblue"),
            (3/25, "blue"),
            (6/25, "lightblue"),
            (9/25, "#ADD8E6"),
            (12/25, "#FFA07A"),
            (15/25, "#FF4500"),
            (18/25, "#FF6347"),
            (21/25, "darkred"),
            (1, "darkred")
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
        kernel_original = GaussianKDE(values_original.T)  # Usando a classe GaussianKDE personalizada
        X_original, Y_original = np.mgrid[self.proj1_data_original.min():self.proj1_data_original.max():100j, 
                                        self.proj2_data_original.min():self.proj2_data_original.max():100j]
        positions_original = np.vstack([X_original.ravel(), Y_original.ravel()])
        Z_original = np.reshape(kernel_original(positions_original.T), X_original.shape)  # Ajuste para compatibilidade
        
        # Adiciona um pequeno valor a Z_original para evitar o erro de divisão por zero
        Z_original = np.clip(Z_original, a_min=1e-10, a_max=None)
        
        G_original = -self.kB * self.temperature * np.log(Z_original)
        G_original = np.clip(G_original - np.min(G_original), 0, 25)
        plt.figure(figsize=(8, 6))
        plt.contourf(X_original, Y_original, G_original, levels=np.linspace(0, 25, 100), cmap=self.custom_cmap)
        plt.colorbar(label='Energia Livre (kJ/mol)', ticks=range(0, 26, 3))
        plt.xlabel('CV1 (Ângulo)')
        plt.ylabel('CV2 (Distância)')
        plt.title('Paisagem Energética Gerada')
        plt.show()


    def main(self):
        self.load_data()
        self.boltzmann_inversion_original(self.proj1_data_original, 'CV1 (Ângulo)')
        self.boltzmann_inversion_original(self.proj2_data_original, 'CV2 (Distância)')
        self.plot_energy_landscape()

if __name__ == "__main__":
    t = 300         # Temperature in K
    kB = 8.314e-3   # Boltzmann constant in kJ/(mol·K)

    cv1_path = './proj1Out.txt'
    cv2_path = './proj2Out.txt'

    fel = FreeEnergyLandscape(cv1_path, cv2_path, t, kB)
    fel.main()
