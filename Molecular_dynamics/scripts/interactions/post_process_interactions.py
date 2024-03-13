import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

class ProcessInteractions:

    def __init__(self, input_dir, output_file):
        self.input_dir = input_dir
        self.output_file = output_file

    def reformat_data(self):
        dat_files = sorted(glob.glob(os.path.join(self.input_dir, 'detailed_interactions_frame_*.dat')))
        with open(self.output_file, 'w') as outfile:
            # Inclui o cabeçalho da coluna de distância
            outfile.write("frame\tligand_atom\tligand_atom_index\taminoacid\taminoacid_atom\taminoacid_number\tdistance\n")
            for dat_file in dat_files:
                with open(dat_file, 'r') as infile:
                    for line in infile:
                        parts = line.strip().split()
                        frame = parts[1].strip(':')
                        ligand_atom = parts[4]
                        ligand_atom_index = parts[6]
                        aminoacid = parts[9]
                        aminoacid_number = parts[10]
                        aminoacid_atom = parts[11]
                        # Extraí a distância do final da linha, assumindo que ela sempre aparece após "Distance:"
                        distance = line.split('Distance: ')[1].strip()
                        # Inclui a distância na linha de saída
                        output_line = f"{frame}\t{ligand_atom}\t{ligand_atom_index}\t{aminoacid}\t{aminoacid_atom}\t{aminoacid_number}\t{distance}\n"
                        outfile.write(output_line)
        print(f"Processamento concluído. Todos os dados foram salvos em: {self.output_file}")



    def calculate_prevalence(self, prevalence_threshold=50):
        data = pd.read_csv(self.output_file, sep='\t')
        print(f"Total de linhas lidas: {data.shape[0]}")  # Diagnóstico: Verificar se os dados estão sendo lidos
        interaction_counts = data.groupby(['ligand_atom', 'ligand_atom_index', 'aminoacid', 'aminoacid_atom', 'aminoacid_number']).size()
        total_frames = data['frame'].nunique()
        print(f"Total de frames únicos: {total_frames}")  # Diagnóstico: Verificar o total de frames
        prevalence = (interaction_counts / total_frames) * 100
        prevalence_filtered = prevalence[prevalence >= prevalence_threshold].reset_index(name='prevalence')
        self.prevalence_data = prevalence_filtered
        if self.prevalence_data.empty:
            print("Nenhuma interação atingiu o limiar de prevalência especificado.")
        else:
            print(self.prevalence_data)


    def plot_prevalence(self):
        if not hasattr(self, 'prevalence_data'):
            print("Prevalence data not calculated. Please run calculate_prevalence first.")
            return

        self.prevalence_data['aminoacid_number'] = pd.to_numeric(self.prevalence_data['aminoacid_number'])
        self.prevalence_data.sort_values(by='aminoacid_number', inplace=True)

        self.prevalence_data['description'] = self.prevalence_data.apply(
            lambda x: f"{x['aminoacid_number']}{x['aminoacid']}({x['aminoacid_atom']}) - {x['ligand_atom']}(LIG)", axis=1)

        plt.figure(figsize=(14, 10))
        plt.bar(self.prevalence_data['description'], self.prevalence_data['prevalence'], color='skyblue')
        plt.ylabel('Prevalência (%)')
        plt.xlabel('Contato Átomo do Aminoácido - Átomo do Ligante')
        plt.xticks(rotation=90, ha="right")
        plt.title('Prevalência de Contato Átomo do Ligante - Átomo do Aminoácido')
        plt.tight_layout()
        plt.show()


    def identify_rotatable_atoms(self, prevalence_threshold=50):
        if not hasattr(self, 'prevalence_data'):
            print("Prevalence data not calculated. Please run calculate_prevalence first.")
            return

        rotatable_atoms = self.prevalence_data[self.prevalence_data['prevalence'] >= prevalence_threshold]
        print(f"Rotatable atoms: \n{rotatable_atoms}")
    
    @staticmethod
    def main():
        input_dir = 'vmd_frames'
        output_file = 'detailed_interaction.tsv'
        prevalence_threshold = 50
        process_interactions = ProcessInteractions(input_dir, output_file)
        process_interactions.reformat_data()
        process_interactions.calculate_prevalence(prevalence_threshold)
        process_interactions.plot_prevalence()
        process_interactions.identify_rotatable_atoms(prevalence_threshold)

if __name__ == '__main__':
    ProcessInteractions.main()