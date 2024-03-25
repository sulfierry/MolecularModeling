from freeEnergyLandscape import *

# Definindo valores padrão
t = 300                     # --temperature           [int] [Kelvin]
kB = 8.314e-3               # --kb                    [float] [kJ/(mol·K)]
energy_threshold = 0.1     # --energy                [float] [kJ/mol]
discrete_val = 0.1         # --discrete              [float]
bins_energy_histogram = 100 # --bins_energy_histogram [int]
kde_bandwidth_cv = None     # --kde_bandwidth         [float]
cv_names = ['CV1', 'CV2']   # --name                  [str] [str]
n_angles = 10               # --gif_angles            [int]
elevation = 10              # --gif_elevation         [int]
duration_per_frame = 0.1    # --gif_duration          [float]
xlim_inf = -8
xlim_sup = 8
ylim_inf = -5
ylim_sup = 7  # Inicialização padrão


class DataProcessor:

    def __init__(self, data_file):
        self.data_file = data_file
        self.index, self.values = self.read_data_file()

    def read_data_file(self):
        index, values = [], []
        with open(self.data_file, 'r') as file:
            for line in file:
                parts = line.strip().split(',')
                value = parts[1].replace('"', '').replace(',', '.')
                index.append(int(parts[0]))
                values.append(float(value))
        return np.array(index), np.array(values)

    def save_processed_data(self, output_file):
        np.savetxt(output_file, np.column_stack((self.index, self.values)), fmt=['%d', '%.6f'], delimiter='\t')

    def plot_energy_distribution(self):
        energy_kde = gaussian_kde(self.values)
        energy_range = np.linspace(self.values.min(), self.values.max(), 500)
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(energy_range, energy_kde(energy_range), color='darkblue', alpha=0.6, lw=2)
        ax.hist(self.values, bins=30, density=True, color='skyblue', alpha=0.5)
        ax.set_title('Data Distribution')
        ax.set_xlabel('Value')
        ax.set_ylabel('Density')
        plt.grid(alpha=0.3)
        plt.show()


    @staticmethod
    def plot_3D_energy_landscape(cv1_path, cv2_path):
        fel = FreeEnergyLandscape(cv1_path, cv2_path, t, kB,
                                  bins=bins_energy_histogram,
                                  kde_bandwidth=kde_bandwidth_cv,
                                  cv_names=cv_names,
                                  discrete=discrete_val,
                                  xlim_inf=xlim_inf, xlim_sup=xlim_sup,
                                  ylim_inf=ylim_inf, ylim_sup=ylim_sup)
        fel.load_data()
        fel.plot_3D_energy_landscape(threshold=energy_threshold, titles=['CV1', 'CV2'])



def main():
    # Processamento e plotagem da distribuição do arquivo de energia
    energy_processor = DataProcessor("./energy_values.csv")
    energy_processor.plot_energy_distribution()
    energy_processor.save_processed_data("./processed_energy_values.txt")

    # Visualizando valores discretos

    # Plotagem da superfície energética em 3D
    DataProcessor.plot_3D_energy_landscape("./processed_energy_values.txt", "./pca_cv1.txt")


if __name__ == "__main__":
    main()
