import os
import re
import subprocess
import matplotlib.pyplot as plt
from collections import defaultdict

class ResiduesPrevalence:

    def __init__(self, ligand_name, distance, threshold=10, color='blue'):
        self.ligand_name = ligand_name
        self.distance = distance
        self.file_path = './residues_near_ligand.dat'
        self.threshold = threshold
        self.color = color

    def search_near_residue(self, topology, trajectory):
        """Chama o VMD para gerar um arquivo com os resíduos próximos a um ligante."""
        tcl_script_filename = 'search_near_residue.tcl'

        tcl_script_content = f'''
        mol load prmtop {topology}.prmtop dcd {trajectory}.dcd
        set ligand_name "{self.ligand_name}"
        set distance {self.distance}
        set closer_residue [atomselect 0 "protein within $distance of resname $ligand_name"]
        set outfile [open "residues_near_ligand.dat" w]
        set num_frames [molinfo 0 get numframes]
        for {{set frame_number 0}} {{$frame_number < $num_frames}} {{incr frame_number}} {{
            $closer_residue frame $frame_number
            $closer_residue update
            set residues [$closer_residue get {{resname resid}}]
            set unique_residues [lsort -unique $residues]
            puts $outfile "Frame $frame_number: [llength $unique_residues] unique residues near $ligand_name"
            foreach {{resname resid}} $unique_residues {{
                puts $outfile "    Residue $resname $resid"
            }}
        }}
        close $outfile
        quit
        '''

        with open(tcl_script_filename, 'w') as file:
            file.write(tcl_script_content)

        subprocess.run(['vmd', '-dispdev', 'text', '-e', tcl_script_filename], check=True)
        os.remove(tcl_script_filename)

    def process_data(self):
        """Processa o arquivo de dados para calcular a prevalência de cada resíduo."""
        residue_count = defaultdict(int)
        with open(self.file_path, 'r') as file:
            for line in file:
                if line.startswith('    Residue '):
                    residues = re.findall(r'([A-Z]+ \d+)', line)
                    for residue in residues:
                        residue_count[residue] += 1

        total_frames = max(int(line.split()[1].rstrip(':')) for line in open(self.file_path) if line.startswith("Frame")) + 1
        self.residue_prevalence = {residue: count / total_frames for residue, count in residue_count.items()}

    def filter_by_threshold(self):
        """Filtra os resíduos com prevalência acima do threshold especificado."""
        self.filtered_data = {res: prevalence for res, prevalence in self.residue_prevalence.items() if prevalence * 100 >= self.threshold}

    def plot_data(self):
        """Gera o gráfico com os dados processados e filtrados, ordenados pelo número do resíduo."""
        sorted_data = sorted(self.filtered_data.items(), key=lambda x: int(x[0].split()[1]))
        labels = [f"{res[0].split()[0]} {res[0].split()[1]}" for res in sorted_data]
        values = [prevalence * 100 for _, prevalence in sorted_data]
        plt.figure(figsize=(10, 8))
        plt.bar(labels, values, color=self.color)
        plt.xlabel('Amino Acid and Number')
        plt.ylabel('Prevalence (%)')
        plt.xticks(rotation=45, ha="right")
        plt.title('Amino Acid Prevalence Near Ligand')
        plt.tight_layout()
        plt.show()

    @staticmethod
    def save_frames_as_pdb(frames, topology, trajectory):
        """Utiliza o VMD para selecionar e salvar uma lista de frames como arquivos PDB."""
        tcl_script_content = f'''
        mol load prmtop {topology}.prmtop dcd {trajectory}.dcd

        foreach frame {{{' '.join(map(str, frames))}}} {{
            animate goto $frame
            set outfile [format "frame%06d.pdb" $frame]
            animate write pdb $outfile
        }}

        quit
        '''

        tcl_script_filename = 'save_frames.tcl'
        with open(tcl_script_filename, 'w') as file:
            file.write(tcl_script_content)

        subprocess.run(['vmd', '-dispdev', 'text', '-e', tcl_script_filename], check=True)
        os.remove(tcl_script_filename)


    @classmethod
    def main(cls):
        # Configurações iniciais: Nome do ligante, distância, topologia e trajetória
        ligand_name = 'LIG'
        distance = 3.0
        topology = 'your_topology'  # Sem a extensão .prmtop
        trajectory = 'your_trajectory'  # Sem a extensão .dcd
        threshold = 10
        color = 'green'

        # Cria uma instância da classe com as configurações desejadas
        instance = cls(ligand_name, distance, threshold, color)

        # Chama o VMD para gerar o arquivo com os resíduos próximos ao ligante
        instance.search_near_residue(topology, trajectory)

        # Processa o arquivo de dados gerado
        instance.process_data()

        # Filtra os dados com base no threshold especificado
        instance.filter_by_threshold()

        # Plota os dados filtrados
        instance.plot_data()

        # Exemplo de frames a serem salvos e chamada para save_frames_as_pdb
        frames_to_save = [0, 10, 20]
        cls.save_frames_as_pdb(frames_to_save, topology, trajectory)


if __name__ == '__main__':

    ResiduesPrevalence.main()
