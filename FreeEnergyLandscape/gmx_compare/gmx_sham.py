import os
import glob
import numpy as np
from PIL import Image
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


class Args:
    def __init__(self):
        self.cv1 = 'proj1Out.txt'
        self.cv2 = 'proj2Out.txt'
        self.o = "GH1-GLUC-AngDist"
        self.t = "Free energy landscape GH1 GLUCOSE"
        self.xl = "Angle"
        self.yl = "Distance"
        self.zl = "Energy (kJ/mol)"
        self.i = "n"
        self.sk = 150
        self.xmin = "30"
        self.xmax = "125"
        self.ymin = "1"
        self.ymax = "20"
        self.septraj = "y"
        self.trj1 = 25000
        self.trj2 = 25000
        self.trj3 = 25000

class FreeEnergyLandscape:
    def __init__(self, args):
        self.args = args

    def read_input_files(self):
        input1 = open(self.args.cv1, 'r').readlines()
        input2 = open(self.args.cv2, 'r').readlines()
        return input1, input2

    def verify_input_sizes(self, input1, input2):
        if len(input1) != len(input2):
            raise ValueError(f'CV1 has different size of CV2 ({len(input1)} x {len(input2)})')

    def generate_table(self, input1, input2):
        table = open(f'{self.args.o}_table_FRAMExCV1xCV2.xvg', 'w+')
        for i1, i2 in zip(input1, input2):
            if not i1.startswith('#') and not i2.startswith('#'):
                frame, cv1 = i1.split()[:2]
                _, cv2 = i2.split()[:2]
                table.write(f"{frame}\t{cv1}\t{cv2}\n")
        table.close()

    def run_gromacs_sham(self):
        # Supondo que {self.args.o}_table_FRAMExCV1xCV2.xvg é o arquivo de entrada correto
        # e que o comando gmx sham irá gerar o arquivo gibbs.xpm
        command = f'gmx sham -f {self.args.o}_table_FRAMExCV1xCV2.xvg -ls gibbs.xpm'
        try:
            os.system(command)
        except Exception as e:
            print(f"Erro ao executar o GROMACS sham: {e}")

    def create_xvg_from_xpm(self, xpm_file):
        out_file = f'{self.args.o}_table_CV1xCV2xENERGY.xvg'
        # Verifica se o arquivo XPM existe antes de tentar abri-lo
        if not os.path.exists(xpm_file):
            raise FileNotFoundError(f"O arquivo {xpm_file} não foi encontrado.")
        
        with open(xpm_file) as xpm_handle:
            xpm_data = []
            x_axis, y_axis, letter_to_value = [], [], {}
            for line in xpm_handle:
                if line.startswith("/* x-axis"):
                    x_axis = [float(x) for x in line.split()[2:-2]]
                elif line.startswith("/* y-axis"):
                    y_axis = [float(y) for y in line.split()[2:-2]]
                elif line.startswith('"') and len(line.split()) > 4:
                    parts = line.split()
                    letter = parts[0][1]
                    value = float(parts[-2].strip('"').split(',')[0])
                    letter_to_value[letter] = value
                elif line.startswith('"') and x_axis and y_axis:
                    xpm_data.insert(0, line.strip().strip(',')[1:-1])

        txt_values = []
        for y_index, row in enumerate(xpm_data):
            for x_index, cell in enumerate(row):
                txt_values.append([x_axis[x_index], y_axis[y_index], letter_to_value[cell]])

        with open(out_file, 'w') as out_handle:
            for x, y, z in txt_values:
                out_handle.write(f"{x}\t{y}\t{z}\n")

    def prepare_data_for_plots(self):
        data = np.genfromtxt(f'{self.args.o}_table_CV1xCV2xENERGY.xvg', delimiter='\t')
        xi, yi, zi = data[:,0], data[:,1], data[:,2]
        xo = np.linspace(xi.min(), xi.max(), 100)
        yo = np.linspace(yi.min(), yi.max(), 100)
        X, Y = np.meshgrid(xo, yo)
        points = np.vstack((xi,yi)).T
        values = zi
        grid_zi = griddata(points, values, (X, Y), method='linear')
        return X, Y, grid_zi, xi, yi

    def make_3d_plot(self, X, Y, Z):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        ax.set_xlabel(self.args.xl)
        ax.set_ylabel(self.args.yl)
        ax.set_zlabel(self.args.zl)
        ax.set_title(self.args.t)
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.savefig(f'{self.args.o}_3Dplot.png')
        plt.close()

    def create_gif(self, X, Y, Z):
        directory = "plot3D_angles"
        if not os.path.exists(directory):
            os.makedirs(directory)
        
        images = []
        for angle in range(0, 360, 10):
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            ax.view_init(30, angle)
            filename = f'{directory}/{self.args.o}_{angle}.png'
            plt.savefig(filename, dpi=96)
            images.append(Image.open(filename))
            plt.close()
        
        gif_filename = f'{self.args.o}_3Dplot.gif'
        images[0].save(gif_filename, save_all=True, append_images=images[1:], duration=100, loop=0)
        # Cleanup
        for filename in glob.glob(f'{directory}/{self.args.o}_*.png'):
            os.remove(filename)

    def make_2d_plot_with_trajectory(self, X, Y, Z, xi, yi):
        plt.figure()
        plt.contourf(X, Y, Z, cmap=cm.coolwarm)
        plt.plot()
        # plt.plot(xi, yi, marker='o', color='black', markersize=2, linestyle="", alpha=0.5)
        plt.colorbar()
        plt.xlabel(self.args.xl)
        plt.ylabel(self.args.yl)
        plt.title(self.args.t)
        plt.savefig(f'{self.args.o}_2Dplot_traj.png')
        plt.close()

    def cleanup_temp_files(self):
        temp_file_patterns = [
            "plot3D_angles/*",
            f"{self.args.o}_table_CV1xCV2xENERGY.xvg",
            f"{self.args.o}_table_FRAMExCV1xCV2.xvg",
            "*.xvg",
            "*.xpm",
            "*.ndx",
        ]
        
        for pattern in temp_file_patterns:
            for file in glob.glob(pattern):
                os.remove(file)
                print(f"Removed: {file}")

    def main(self):
        
        input1, input2 = self.read_input_files()
        self.verify_input_sizes(input1, input2)
        self.generate_table(input1, input2)
        self.run_gromacs_sham()
        self.create_xvg_from_xpm('gibbs.xpm')
        X, Y, Z, xi, yi = self.prepare_data_for_plots()
        self.make_3d_plot(X, Y, Z)
        self.create_gif(X, Y, Z)
        if self.args.septraj == 'y':
            self.make_2d_plot_with_trajectory(X, Y, Z, xi, yi)
        self.cleanup_temp_files()

if __name__ == "__main__":
    args = Args() 
    fel = FreeEnergyLandscape(args)
    fel.main()

