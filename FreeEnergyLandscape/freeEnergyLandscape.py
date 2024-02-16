import os
import sys
import shutil
import tempfile
import numpy as np
import imageio.v2 as imageio
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
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


    def boltzmann_inversion_original(self, data, title, threshold=10):
        # Calcular G e plotar o gráfico de energia livre
        hist, bin_edges = np.histogram(data, bins=100, density=True)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        hist = np.clip(hist, a_min=1e-10, a_max=None)
        G = -self.kB * self.temperature * np.log(hist / np.max(hist))

        plt.figure(figsize=(10, 6))
        plt.plot(bin_centers, G, label='Free energy', color='red')

        # Para plotar as bolinhas magentas, precisamos identificar os bins de baixa energia
        low_energy_bins = G < threshold
        if any(low_energy_bins):
            for bin_center, g_value in zip(bin_centers[low_energy_bins], G[low_energy_bins]):
                # Achar os índices dos dados que caem dentro deste bin
                indices = np.where((data >= bin_center - (bin_centers[1]-bin_centers[0])/2) & 
                                (data < bin_center + (bin_centers[1]-bin_centers[0])/2))[0]
                # Plotar as bolinhas magentas para cada valor de CV neste bin
                plt.scatter([bin_center] * len(indices), [g_value] * len(indices), color='magenta')

        plt.xlabel(title)
        plt.ylabel('Free energy (kJ/mol)')
        plt.grid(True)
        plt.legend()
        plt.title(f'Free Energy Profile for {title}')
        plt.show()
        
        self.identify_and_save_low_energy_frames(data, G, bin_centers, threshold, title)


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
        """
        Calcula a energia livre e prepara os dados para plotagem.
        :param data: Dados de entrada para os quais a energia livre será calculada.
        :return: Dicionário contendo 'X_original', 'Y_original', e 'G_original' para plotagem.
        """
        values_original = np.vstack([data[:, 0], data[:, 1]]).T
        kernel_original = gaussian_kde(values_original.T)
        X_original, Y_original = np.mgrid[data[:, 0].min():data[:, 0].max():100j, 
                                        data[:, 1].min():data[:, 1].max():100j]
        positions_original = np.vstack([X_original.ravel(), Y_original.ravel()])
        Z_original = np.reshape(kernel_original(positions_original).T, X_original.shape)
        G_original = -self.kB * self.temperature * np.log(Z_original)
        G_original = np.clip(G_original - np.min(G_original), 0, 25)
        
        return {'X_original': X_original, 'Y_original': Y_original, 'G_original': G_original}


    def plot_energy_landscape(self, threshold=10):
        # Combina os dados das variáveis coletivas
        data = np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None]))
        result = self.calculate_free_energy(data)
        plt.figure(figsize=(8, 6))

        # Plotagem da paisagem de energia
        cont = plt.contourf(result['X_original'], result['Y_original'], result['G_original'], 
                            levels=np.linspace(np.min(result['G_original']), np.max(result['G_original']), 100), 
                            cmap=self.custom_cmap, extend='both')

        # Identificação dos pontos de baixa energia
        low_energy_mask = result['G_original'] <= threshold
        plt.scatter(result['X_original'][low_energy_mask], result['Y_original'][low_energy_mask], 
                    color='magenta', s=10, label=f'Low Energy (<= {threshold} kJ/mol)')

        plt.legend(loc='lower left', bbox_to_anchor=(1, 1))
        cbar = plt.colorbar(cont)
        cbar.set_label('Free energy (kJ/mol)')
        plt.xlabel('CV1 (Angle)')
        plt.ylabel('CV2 (Distance)')
        plt.title('Free energy landscape')
        plt.show()

        self.save_low_energy_frames_to_tsv_2D(threshold, 'low_energy_frames.tsv')

        # Salvar os pontos de baixa energia e os frames correspondentes
        with open('landscape_cv1_cv2_minimuns.tsv', 'w') as file:
            file.write("frame\tcv1\tcv2\tenergy\n")
            # Assumindo que cada ponto em data corresponde a um frame, começando do 0
            for i in np.where(low_energy_mask.ravel())[0]:
                frame = i  # Diretamente associado, pois cada linha em data corresponde a um frame
                energy = result['G_original'].ravel()[i]
                cv1 = result['X_original'].ravel()[i]
                cv2 = result['Y_original'].ravel()[i]
                file.write(f"{frame}\t{cv1}\t{cv2}\t{energy}\n")


        # Salvamento dos pontos de baixa energia
        low_energy_x, low_energy_y = result['X_original'][low_energy_mask], result['Y_original'][low_energy_mask]
        low_energy_g = result['G_original'][low_energy_mask]
        np.savetxt('landscape_cv1_cv2_minimuns.tsv', np.column_stack((low_energy_x.ravel(), low_energy_y.ravel(), low_energy_g.ravel())), fmt='%f', header='CV1\tCV2\tEnergy', comments='')


    def save_low_energy_points_to_tsv(self, threshold=10):
        # Calcula a paisagem de energia
        data = np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None]))
        result = self.calculate_free_energy(data)

        # Identifica os pontos de baixa energia (<= threshold)
        low_energy_mask = result['G_original'] <= threshold

        # Prepara os dados para salvar: coordenadas X e Y dos pontos de baixa energia e seus valores de energia
        low_energy_x, low_energy_y = np.meshgrid(np.linspace(data[:, 0].min(), data[:, 0].max(), 100), 
                                                np.linspace(data[:, 1].min(), data[:, 1].max(), 100))
        low_energy_g = result['G_original'][low_energy_mask]

        # Achata os arrays para simplificar o uso
        low_energy_x_flat = low_energy_x[low_energy_mask].ravel()
        low_energy_y_flat = low_energy_y[low_energy_mask].ravel()
        low_energy_g_flat = low_energy_g.ravel()

        # Abre o arquivo para escrita
        with open('landscape_cv1_cv2_minimuns.tsv', 'w') as file:
            file.write("cv1\tcv2\tenergy\n")
            
            # Loop sobre os pontos de baixa energia
            for x, y, g in zip(low_energy_x_flat, low_energy_y_flat, low_energy_g_flat):
                file.write(f"{x}\t{y}\t{g}\n")
        
        print("Low energy points saved to landscape_cv1_cv2_minimuns.tsv")

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



    def main(self):
        self.load_data()
        self.boltzmann_inversion_original(self.proj1_data_original, 'CV1 (Angle)')
        self.boltzmann_inversion_original(self.proj2_data_original, 'CV2 (Distance)')
        self.plot_energy_landscape()
        self.save_low_energy_points_to_tsv(threshold=5)
        # self.plot_3D_energy_landscape()
        # self.create_3D_gif()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python freeEnergyLandscape.py path/to/cv1_data.txt path/to/cv2_data.txt")
        sys.exit(1)

    cv1_path, cv2_path = sys.argv[1:3]

    try:
        t = 300  # Temperature in K
        kB = 8.314e-3  # Boltzmann constant in kJ/(mol·K)
        
        fel = FreeEnergyLandscape(cv1_path, cv2_path, t, kB)
        fel.main()


    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

