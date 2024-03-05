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
        plt.savefig('aac_prevalence.png')
        plt.tight_layout()
        plt.show()



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

        # Processa o arquivo de dados gerado
        instance.process_data()

        # Filtra os dados com base no threshold especificado
        instance.filter_by_threshold()

        # Plota os dados filtrados
        instance.plot_data()


if __name__ == '__main__':

    ResiduesPrevalence.main()
