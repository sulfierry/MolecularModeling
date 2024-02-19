import os
import sys
import shutil
import argparse
import tempfile
import numpy as np
import imageio.v2 as imageio
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from matplotlib.colors import LinearSegmentedColormap

class FreeEnergyLandscape:
    
    def __init__(self, cv1_path, cv2_path, temperature, boltzmann_constant, bins=100, kde_bandwidth=None):
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
        self.bins = bins
        self.kde_bandwidth = kde_bandwidth

    def load_data(self):
        self.proj1_data_original = np.loadtxt(self.cv1_path, usecols=[1])
        self.proj2_data_original = np.loadtxt(self.cv2_path, usecols=[1])

    def boltzmann_inversion(self, data_list, titles, threshold):
        fig_combined, axs_combined = plt.subplots(1, len(data_list), figsize=(20, 6), sharey=True)

        # Primeiro, calculamos o G_max_global com base em todos os conjuntos de dados
        G_max_global = -np.inf
        for data in data_list:
            hist, bin_edges = np.histogram(data, bins=100, density=True)
            hist = np.clip(hist, a_min=1e-10, a_max=None)
            G = -self.kB * self.temperature * np.log(hist + 1e-10)
            G_max_global = max(G_max_global, np.max(G))

        # Agora, plotamos cada variável coletiva e seus pontos de baixa energia
        for ax, data, title in zip(axs_combined, data_list, titles):
            hist, bin_edges = np.histogram(data, bins=100, density=True)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            G = -self.kB * self.temperature * np.log(hist + 1e-10)
            G_normalized = G / G_max_global

            ax.plot(bin_centers, G_normalized, label='Free energy', color='red')
            ax.set_xlabel(title)
            ax.set_ylim(0, 1)

            # Aplicamos o threshold diretamente ao G normalizado para cada variável
            if threshold is not None:
                # Para cada variável, recalculamos o threshold_normalized baseado no seu próprio G_max
                G_max = np.max(G)
                threshold_normalized = threshold / G_max
                low_energy_bins = G_normalized < threshold_normalized
                if any(low_energy_bins):
                    ax.scatter(bin_centers[low_energy_bins], G_normalized[low_energy_bins], color='magenta', label='Low energy points')

        axs_combined[0].set_ylabel('Normalized Free Energy')
        plt.legend()
        plt.suptitle('Normalized Free Energy Profile Comparison')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.savefig('Combined_Free_Energy_Profile_Normalized.png')
        plt.show()


    def identify_and_save_low_energy_frames(self, data, G, bin_centers, threshold, title):
        hist, bin_edges = np.histogram(data, bins=100, density=True)
        bin_width = bin_edges[1] - bin_edges[0]
        filename = f"{title.replace(' ', '_').lower()}_low_energy_frames.tsv"

        with open(filename, 'w') as file:
            file.write("frame\tcv\tenergy\n")
            for bin_center, energy in zip(bin_centers, G):
                if energy < threshold:
                    # Identifica índices dos dados que caem dentro deste bin de baixa energia
                    indices_within_bin = np.where((data >= bin_center - bin_width / 2) & 
                                                (data < bin_center + bin_width / 2))[0]
                    
                    for frame in indices_within_bin:
                        file.write(f"{frame}\t{data[frame]}\t{energy}\n")
        print(f"Saved low energy frames to {filename}")



    def calculate_free_energy(self, data):
        if hasattr(self, 'cached_results'):
            return self.cached_results

        values_original = np.vstack([data[:, 0], data[:, 1]]).T
        if self.kde_bandwidth:
            kernel_original = gaussian_kde(values_original.T, bw_method=self.kde_bandwidth)
        else:
            kernel_original = gaussian_kde(values_original.T)
        X_original, Y_original = np.mgrid[data[:, 0].min():data[:, 0].max():100j, 
                                          data[:, 1].min():data[:, 1].max():100j]
        positions_original = np.vstack([X_original.ravel(), Y_original.ravel()])
        Z_original = np.reshape(kernel_original(positions_original).T, X_original.shape)
        G_original = -self.kB * self.temperature * np.log(Z_original)
        G_original = np.clip(G_original - np.min(G_original), 0, 25)

        self.cached_results = {'X_original': X_original, 'Y_original': Y_original, 'G_original': G_original}
        return self.cached_results


    def plot_energy_landscape(self, threshold):
        data = np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None]))
        result = self.calculate_free_energy(data)
        plt.figure(figsize=(8, 6))

        cont = plt.contourf(result['X_original'], result['Y_original'], result['G_original'], 
                            levels=np.linspace(np.min(result['G_original']), np.max(result['G_original']), 100), 
                            cmap=self.custom_cmap, extend='both')

        if threshold is not None:
            low_energy_mask = result['G_original'] <= threshold
            if np.any(low_energy_mask):  # Verifica se existem pontos de baixa energia para plotar
                plt.scatter(result['X_original'][low_energy_mask], result['Y_original'][low_energy_mask], 
                            color='magenta', s=10, label=f'Energy <= {threshold} kJ/mol')
                plt.legend(loc='lower left', bbox_to_anchor=(1, 1))

        # Salvar os pontos de baixa energia e os frames correspondentes sempre, independente da visualização
        if threshold is not None:
            with open('landscape_cv1_cv2_minimuns.tsv', 'w') as file:
                file.write("frame\tcv1\tcv2\tenergy\n")
                for i in np.where(low_energy_mask.ravel())[0]:
                    # Nota: Esta parte assume uma correspondência direta que pode precisar de revisão para precisão
                    energy = result['G_original'].ravel()[i]
                    cv1 = result['X_original'].ravel()[i]
                    cv2 = result['Y_original'].ravel()[i]
                    # A identificação precisa do 'frame' corresponde a um desafio adicional aqui
                    file.write(f"{i}\t{cv1}\t{cv2}\t{energy}\n")

        cbar = plt.colorbar(cont)
        cbar.set_label('Free energy (kJ/mol)')
        plt.xlabel('CV1 (Angle)')
        plt.ylabel('CV2 (Distance)')
        plt.title('Free energy landscape')
        plt.show()

    def save_low_energy_frames_to_tsv_2D(self, threshold=5, filename='low_energy_frames.tsv'):
            data = np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None]))
            result = self.calculate_free_energy(data)

            # Mascara para identificar pontos de energia abaixo do limiar
            low_energy_mask = result['G_original'] <= threshold

            # Preparar e salvar os dados
            with open(filename, 'w') as file:
                file.write("frame\tcv1\tcv2\tenergy\n")
                for idx in np.where(low_energy_mask.ravel())[0]:
                    # Encontra o frame mais próximo para cada ponto de baixa energia
                    x, y = np.meshgrid(np.linspace(data[:, 0].min(), data[:, 0].max(), 100), 
                                    np.linspace(data[:, 1].min(), data[:, 1].max(), 100), indexing='ij')
                    cv1_val, cv2_val = x.ravel()[idx], y.ravel()[idx]
                    energy_val = result['G_original'].ravel()[idx]

                    # Identificação aproximada do frame
                    frame = np.argmin((self.proj1_data_original - cv1_val)**2 + (self.proj2_data_original - cv2_val)**2)
                    file.write(f"{frame}\t{cv1_val}\t{cv2_val}\t{energy_val}\n")


    def save_low_energy_points_to_tsv(self, threshold):
        if threshold is None or threshold < 0:
            print("Threshold is None or less than 0. No points will be saved.")
            return

        # Carrega os dados originais novamente para obter os números dos frames
        cv1_data = np.loadtxt(self.cv1_path)
        cv2_data = np.loadtxt(self.cv2_path)

        # Supondo que a primeira coluna de cada arquivo contém os números dos frames
        frames_cv1 = cv1_data[:, 0].astype(int)
        frames_cv2 = cv2_data[:, 0].astype(int)

        # Calcula a paisagem de energia usando os dados carregados previamente
        result = self.calculate_free_energy(np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None])))

        low_energy_mask = result['G_original'] <= threshold
        if not np.any(low_energy_mask):
            print("No low energy points found. No points will be saved.")
            return

        with open('low_energy_points.tsv', 'w') as file:
            file.write("frame\tcv1\tcv2\tenergy\n")

            # Itera sobre os pontos de baixa energia
            for idx in np.where(low_energy_mask.ravel())[0]:
                # A coordenada X, Y no espaço de CV1 e CV2 pode ser mapeada de volta para o frame mais próximo
                # Aqui, assumimos que os frames em CV1 e CV2 estão alinhados e usamos frames de CV1 como referência
                frame = frames_cv1[idx]  # Usa o número do frame diretamente dos dados carregados
                cv1 = self.proj1_data_original[idx]
                cv2 = self.proj2_data_original[idx]
                g = result['G_original'].flatten()[idx]
                file.write(f"{frame}\t{cv1}\t{cv2}\t{g}\n")

        print("Low energy points saved to low_energy_points.tsv")



    def plot_3D_energy_landscape(self):
        data = np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None]))
        result = self.calculate_free_energy(data)
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(result['X_original'], result['Y_original'], result['G_original'], cmap=self.custom_cmap, edgecolor='none')
        ax.set_xlabel('CV1 (Angle)')
        ax.set_ylabel('CV2 (Distance)')
        ax.set_title('3D Free energy landscape')
        fig.colorbar(surf, shrink=0.5, aspect=5, label='Free energy (kJ/mol)')
        plt.show()

    def create_3D_gif(self, gif_filename='energy_landscape_3D.gif', n_angles=10, elevation=15, duration_per_frame=0.01):
        temp_dir = tempfile.mkdtemp()  # Cria um diretório temporário para armazenar os frames
        filenames = []

        # Utiliza a função calculate_free_energy refatorada para obter os dados
        # Aqui, combinamos proj1_data_original e proj2_data_original para a visualização
        data = np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None]))
        result = self.calculate_free_energy(data)

        # Gera uma lista de ângulos para um movimento contínuo e suave
        angles = list(range(0, 360, int(360 / n_angles)))

        for i, angle in enumerate(angles):
            fig = plt.figure(figsize=(10, 7))
            ax = fig.add_subplot(111, projection='3d')
            surf = ax.plot_surface(result['X_original'], result['Y_original'], result['G_original'], cmap=self.custom_cmap, edgecolor='none', vmin=np.min(result['G_original']), vmax=np.max(result['G_original']))
            ax.view_init(elev=elevation, azim=angle)
            ax.set_xlabel('CV1 (Angle)')
            ax.set_ylabel('CV2 (Distance)')
            ax.set_zlabel('Free energy (kJ/mol)')
            ax.set_title('3D Free energy landscape')

            # Adiciona a barra de cores no primeiro e último frame
            if i == 0 or i == len(angles) - 1:
                fig.colorbar(surf, shrink=0.5, aspect=5, label='Free energy (kJ/mol)')

            frame_filename = os.path.join(temp_dir, f"frame_{angle:03d}.png")
            plt.savefig(frame_filename)
            filenames.append(frame_filename)
            plt.close()

        with imageio.get_writer(gif_filename, mode='I', duration=duration_per_frame) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)

        shutil.rmtree(temp_dir)  # Limpa os arquivos temporários

    def plot_histogram(self, data_list, titles):
        # Plotando histogramas absolutos individualmente
        for data, title in zip(data_list, titles):
            plt.figure(figsize=(8, 6))
            plt.hist(data, bins=self.bins, density=True, alpha=0.7, color='blue')
            plt.title(f'Absolute {title}')
            plt.xlabel('Value')
            plt.ylabel('Frequency')
            plt.grid(True)
            plt.savefig(f'histogram_absolute_{title.replace(" ", "_")}.png')
            plt.close()  # Fecha a figura após salvar

        # Plotando histogramas normalizados lado a lado
        plt.figure(figsize=(8 * len(data_list), 6))
        for i, (data, title) in enumerate(zip(data_list, titles)):
            # Normalização dos dados
            data_normalized = (data - np.min(data)) / (np.max(data) - np.min(data))
            
            # Criação de subplots lado a lado
            ax = plt.subplot(1, len(data_list), i + 1)
            ax.hist(data_normalized, bins=self.bins, density=True, alpha=0.7, color='green')
            ax.set_title(f'Normalized {title}')
            ax.set_xlabel('Normalized Value')
            if i == 0:  # Apenas o primeiro subplot terá o label do eixo Y
                ax.set_ylabel('Frequency')
            ax.grid(True)

        plt.tight_layout()
        plt.savefig('histograms_normalized_side_by_side.png')
        plt.show()


    def cv_by_frame(self, data_list, titles):
        frames = np.arange(len(data_list[0]))  # Assumindo que todos os conjuntos têm o mesmo número de frames
        for data, title in zip(data_list, titles):
            plt.figure(figsize=(10, 6))
            plt.plot(frames, data, label=title)
            plt.xlabel('Frame')
            plt.ylabel(title)
            plt.title(f'CV by Frame - {title}')
            plt.legend()
            plt.savefig(f'cv_by_frame_absolute_{title.replace(" ", "_")}.png')
            plt.close()

        # Plot combinado dos valores relativos
        plt.figure(figsize=(10, 6))
        for data, title in zip(data_list, titles):
            data_normalized = (data - np.min(data)) / (np.max(data) - np.min(data))
            plt.plot(frames, data_normalized, label=title)
        plt.xlabel('Frame')
        plt.ylabel('Normalized CV')
        plt.title('CV by Frame - Combined Normalized')
        plt.legend()
        plt.savefig('cv_by_frame_combined_normalized.png')
        plt.show()


    def verify_input(self, data_path):
        try:
            data = np.loadtxt(data_path, dtype=float)
            if data.shape[1] != 2:
                raise ValueError("O arquivo de entrada deve ter exatamente duas colunas.")

            frames = data[:, 0]
            if not np.all(frames == np.arange(len(frames))):
                raise ValueError("A primeira coluna deve conter valores inteiros sequenciais que representam os frames.")

        except ValueError as e:
            print(f"Erro na verificação do arquivo {data_path}: {e}")
            sys.exit(1)


    def main(self, energy_threshold):
        
        # Verificar ambos os arquivos de entrada antes de carregar os dados
        self.verify_input(self.cv1_path)
        self.verify_input(self.cv2_path)
        self.load_data()
        
        self.boltzmann_inversion(
            data_list=[self.proj1_data_original, self.proj2_data_original], 
            titles=['CV1 (Angle)', 'CV2 (Distance)'], 
            threshold=energy_threshold)
        
        self.plot_histogram(
            data_list=[self.proj1_data_original, self.proj2_data_original], 
            titles=['CV1 (Angle)', 'CV2 (Distance)'])

        self.cv_by_frame(
            data_list=[self.proj1_data_original, self.proj2_data_original], 
            titles=['CV1 (Angle)', 'CV2 (Distance)'])


        self.plot_energy_landscape(threshold=energy_threshold)

        # self.plot_3D_energy_landscape()
        # self.create_3D_gif()
        
        self.save_low_energy_points_to_tsv(threshold=energy_threshold) # save low energy frames to tsv

    
        # Após o uso final dos dados, limpe-os para liberar memória
        if hasattr(self, 'cached_results'):
            del self.cached_results



if __name__ == "__main__":

    # Definindo valores padrão
    t = 300                      # --temperature           [int] [Kelvin]
    kB = 8.314e-3                # --kb                    [float] [kJ/(mol·K)]
    energy_threshold = None      # --energy                [int] [kJ/mol]
    bins_energy_histogram = 100  # --bins_energy_histogram [int]
    kde_bandwidth_cv = None      # --kde_bandwidth         [float]


    if len(sys.argv) >= 3:
        cv1_path, cv2_path = sys.argv[1], sys.argv[2]

        # Processar argumentos adicionais como pares chave-valor
        for i in range(3, len(sys.argv), 2):
            if i+1 < len(sys.argv):
                key = sys.argv[i]
                value = sys.argv[i+1]
                if key == "--temperature":
                    t = float(value)
                elif key == "--kb":
                    kB = float(value)
                elif key == "--energy":
                    energy_threshold = float(value)
                elif key == "--bins_energy_histogram":
                    bins_energy_histogram = int(value)
                elif key == "--kde_bandwidth":
                    if value.lower() != "none":
                        kde_bandwidth_cv = float(value)
                    else:
                        kde_bandwidth_cv = None
    else:
        print("Usage: python freeEnergyLandscape.py path/to/cv1_data.txt path/to/cv2_data.txt [optional arguments --temperature, --kb, --energy_threshold, --bins_energy_histogram, --kde_bandwidth_cv]")
        sys.exit(1)

    try:
        fel = FreeEnergyLandscape(cv1_path, cv2_path, t, kB,  
                                bins=bins_energy_histogram, 
                                kde_bandwidth=kde_bandwidth_cv)
        fel.main(energy_threshold)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)
