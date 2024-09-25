import os
import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt

class ProcessInteractions:

    def __init__(self, input_dir, output_file):
        self.input_dir = input_dir
        self.output_file = output_file
        self.data = None

    def load_data(self):
        dat_files = sorted(glob.glob(os.path.join(self.input_dir, 'detailed_interactions_frame_*.dat')))
        print(f"Número de arquivos .dat encontrados: {len(dat_files)}")
        
        data_list = []
        for dat_file in dat_files:
            with open(dat_file, 'r') as infile:
                for line in infile:
                    parts = line.strip().split()
                    frame = int(parts[1].strip(':'))
                    data_list.append({
                        'frame': frame,
                        'chainA_residue': parts[4],
                        'chainA_residue_number': int(parts[5]),
                        'chainA_atom': parts[6],
                        'chainB_residue': parts[11],
                        'chainB_residue_number': int(parts[12]),
                        'chainB_atom': parts[13],
                        'distance': float(parts[-1])
                    })
        
        self.data = pd.DataFrame(data_list)
        print(f"Total de interações carregadas: {len(self.data)}")

    def calculate_prevalence(self, prevalence_threshold):
        if self.data is None:
            print("Dados não carregados. Execute load_data() primeiro.")
            return

        print(f"Dimensões do DataFrame: {self.data.shape}")
        print("Primeiras linhas dos dados:")
        print(self.data.head())
        
        total_frames = self.data['frame'].nunique()
        print(f"\nNúmero total de frames: {total_frames}")
        
        interaction_counts = self.data.groupby(['chainA_residue', 'chainA_residue_number', 'chainA_atom', 'chainB_residue', 'chainB_residue_number', 'chainB_atom']).size()
        prevalence = (interaction_counts / total_frames) * 100
        
        prevalence_filtered = prevalence[prevalence >= prevalence_threshold].reset_index(name='prevalence')
        self.prevalence_data = prevalence_filtered
        
        print(f"\nNúmero de interações que atingiram o limiar de prevalência ({prevalence_threshold}%): {len(self.prevalence_data)}")
        
        if self.prevalence_data.empty:
            print("Nenhuma interação atingiu o limiar de prevalência especificado.")
        else:
            print("\nInterações que atingiram o limiar de prevalência:")
            print(self.prevalence_data)

        # Filtrar os dados originais com base nas interações que atingiram o limiar de prevalência
        self.filtered_data = self.data.merge(self.prevalence_data, on=['chainA_residue', 'chainA_residue_number', 'chainA_atom', 'chainB_residue', 'chainB_residue_number', 'chainB_atom'])
        
        # Salvar os dados filtrados
        self.filtered_data.to_csv(self.output_file, sep='\t', index=False)
        print(f"\nDados filtrados salvos em: {self.output_file}")
        print(f"Número de linhas no arquivo de saída: {len(self.filtered_data)}")

    def plot_prevalence(self):
        if not hasattr(self, 'prevalence_data'):
            print("Dados de prevalência não calculados. Por favor, execute calculate_prevalence primeiro.")
            return
        
        self.prevalence_data['chainB_residue_number'] = pd.to_numeric(self.prevalence_data['chainB_residue_number'], errors='coerce')
        self.prevalence_data = self.prevalence_data.sort_values(by='chainB_residue_number', na_position='last')
        
        self.prevalence_data['description'] = self.prevalence_data.apply(
            lambda x: f"{x['chainB_residue_number']}{x['chainB_residue']}({x['chainB_atom']}) - {x['chainA_residue_number']}{x['chainA_residue']}({x['chainA_atom']})", axis=1)
        
        plt.figure(figsize=(14, 10))
        plt.bar(self.prevalence_data['description'], self.prevalence_data['prevalence'], color='skyblue')
        plt.ylabel('Prevalência (%)')
        plt.xlabel('Contato Átomo da Cadeia B - Átomo da Cadeia A')
        plt.xticks(rotation=90, ha="right")
        plt.title('Prevalência de Contato Átomo da Cadeia A - Átomo da Cadeia B')
        plt.tight_layout()
        plt.savefig('prevalence.png')
        print("Gráfico de prevalência salvo como 'prevalence.png'")
        plt.show()

    def identify_rotatable_atoms(self, prevalence_threshold):
        if not hasattr(self, 'prevalence_data'):
            print("Dados de prevalência não calculados. Por favor, execute calculate_prevalence primeiro.")
            return
        rotatable_atoms = self.prevalence_data[self.prevalence_data['prevalence'] >= prevalence_threshold]
        print(f"Átomos rotacionáveis (prevalência >= {prevalence_threshold}%):")
        print(rotatable_atoms)
        print(f"Total de átomos rotacionáveis identificados: {len(rotatable_atoms)}")

if __name__ == '__main__':
    input_dir = 'vmd_frames'
    output_file = 'detailed_interaction50.tsv'
    prevalence_threshold = 50  # default value
    if len(sys.argv) > 1:
        prevalence_threshold = float(sys.argv[1])
    
    print(f"Diretório de entrada: {input_dir}")
    print(f"Arquivo de saída: {output_file}")
    print(f"Limiar de prevalência: {prevalence_threshold}%")
    
    process_interactions = ProcessInteractions(input_dir, output_file)
    process_interactions.load_data()
    process_interactions.calculate_prevalence(prevalence_threshold)
    #process_interactions.plot_prevalence()
    process_interactions.identify_rotatable_atoms(prevalence_threshold)
