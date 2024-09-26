# apenas tenha a certeza que a pasta 'vmd_frames' esta no diretorio atual
# esta pasta é obtida apos a execucao do script 'near_atom_from_chain.sh'

import os
import sys
import glob
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


class ProcessInteractions:

    def __init__(self, input_dir, output_file, prevalence_threshold):
        self.input_dir = input_dir
        self.output_file = output_file
        self.prevalence_threshold = prevalence_threshold
        self.data = None
        self.prevalence_data = None
        self.filtered_data = None

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
        print("Átomos únicos da cadeia B nos dados carregados:")
        print(self.data['chainB_atom'].unique())

    def calculate_prevalence(self):
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
        
        self.prevalence_data = prevalence[prevalence >= self.prevalence_threshold].reset_index(name='prevalence')
        
        print(f"\nNúmero de interações que atingiram o limiar de prevalência ({self.prevalence_threshold}%): {len(self.prevalence_data)}")
        
        if self.prevalence_data.empty:
            print("Nenhuma interação atingiu o limiar de prevalência especificado.")
        else:
            print("\nInterações que atingiram o limiar de prevalência:")
            print(self.prevalence_data)

        print("\nÁtomos únicos da cadeia B nos dados de prevalência:")
        print(self.prevalence_data['chainB_atom'].unique())

        # Filtrar os dados originais com base nas interações que atingiram o limiar de prevalência
        self.filtered_data = self.data.merge(self.prevalence_data, on=['chainA_residue', 'chainA_residue_number', 'chainA_atom', 'chainB_residue', 'chainB_residue_number', 'chainB_atom'])
        
        print("\nÁtomos únicos da cadeia B nos dados filtrados:")
        print(self.filtered_data['chainB_atom'].unique())

        # Verificação explícita dos átomos da cadeia B
        unique_chainB_atoms = self.filtered_data['chainB_atom'].unique()
        print(f"Átomos únicos da cadeia B nos dados filtrados: {unique_chainB_atoms}")
        if len(unique_chainB_atoms) == 1 and unique_chainB_atoms[0] == 'N':
            print("AVISO: Apenas o átomo 'N' está presente para a cadeia B nos dados filtrados.")

        # Salvar os dados filtrados
        self.filtered_data.to_csv(self.output_file, sep='\t', index=False)
        print(f"\nDados filtrados salvos em: {self.output_file}")
        print(f"Número de linhas no arquivo de saída: {len(self.filtered_data)}")

    def plot_prevalence(self):
        if self.prevalence_data is None:
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

    def identify_rotatable_atoms(self):
        if self.prevalence_data is None:
            print("Dados de prevalência não calculados. Por favor, execute calculate_prevalence primeiro.")
            return
        rotatable_atoms = self.prevalence_data[self.prevalence_data['prevalence'] >= self.prevalence_threshold]
        print(f"Átomos rotacionáveis (prevalência >= {self.prevalence_threshold}%):")
        print(rotatable_atoms)
        print(f"Total de átomos rotacionáveis identificados: {len(rotatable_atoms)}")

    def run(self):
        self.load_data()
        self.calculate_prevalence()
        self.identify_rotatable_atoms()
        # Uncomment the following line if you want to plot the prevalence
        # self.plot_prevalence()

    @staticmethod
    def main():
        parser = argparse.ArgumentParser(description='Process protein chain interactions.')
        parser.add_argument('--input_dir', type=str, default='vmd_frames', help='Directory containing input .dat files')
        parser.add_argument('--output_file', type=str, default='detailed_interactions.tsv', help='Output file name')
        parser.add_argument('--prevalence_threshold', type=float, default=50, help='Prevalence threshold for interactions')
        args = parser.parse_args()

        print(f"Diretório de entrada: {args.input_dir}")
        print(f"Arquivo de saída: {args.output_file}")
        print(f"Limiar de prevalência: {args.prevalence_threshold}%")

        process_interactions = ProcessInteractions(args.input_dir, args.output_file, args.prevalence_threshold)
        process_interactions.run()

    
class ViewInteractions:
    def __init__(self, input_file, top_interactions):
        self.input_file = input_file
        self.top_interactions = top_interactions
        self.data = None
        self.interaction_grouped = None

    def load_data(self):
        if not os.path.exists(self.input_file):
            print(f"Erro: O arquivo {self.input_file} não existe.")
            exit(1)
        try:
            self.data = pd.read_csv(self.input_file, sep='\t')
        except Exception as e:
            print(f"Erro ao ler o arquivo: {e}")
            exit(1)

    def process_data(self):
        print("Primeiras linhas do arquivo de entrada:")
        print(self.data.head())

        print("\nÁtomos únicos na cadeia A (antes do agrupamento):")
        print(self.data['chainA_atom'].unique())
        print("\nÁtomos únicos na cadeia B (antes do agrupamento):")
        print(self.data['chainB_atom'].unique())

        print("\nContagem de átomos na cadeia B:")
        print(self.data['chainB_atom'].value_counts())

        self.data['interaction'] = (self.data['chainA_residue'] + ' (' + self.data['chainA_residue_number'].astype(str) + ') - ' +
                                    self.data['chainB_residue'] + ' (' + self.data['chainB_residue_number'].astype(str) + ')')

        self.interaction_grouped = self.data.groupby('interaction').apply(self.aggregate_interaction_data).reset_index()
        self.interaction_grouped = self.interaction_grouped.sort_values(by='prevalence', ascending=False)

        print("\nPrimeiras interações mais prevalentes:")
        print(self.interaction_grouped.head())

        print("\nÁtomos únicos na cadeia A (após o agrupamento):")
        print(self.interaction_grouped['chainA_atom'].unique())
        print("\nÁtomos únicos na cadeia B (após o agrupamento):")
        print(self.interaction_grouped['chainB_atom'].unique())

        print("\nAmostra de linhas para chainB_atom:")
        print(self.interaction_grouped[['interaction', 'chainB_atom']].sample(min(10, len(self.interaction_grouped))))

    @staticmethod
    def aggregate_interaction_data(group):
        chainB_atoms = group['chainB_atom'].value_counts()
        most_common_chainB_atom = chainB_atoms.index[0] if not chainB_atoms.empty else 'Unknown'
        return pd.Series({
            'prevalence': group['prevalence'].mean(),
            'distance': group['distance'].mean(),
            'chainA_atom': group['chainA_atom'].mode().iloc[0],
            'chainB_atom': most_common_chainB_atom
        })

    def save_csv(self):
        output_csv = 'aminoacid_interactions_prevalence.csv'
        self.interaction_grouped.to_csv(output_csv, index=False)
        print(f"\nArquivo CSV salvo como '{output_csv}'")

    def plot_interactions(self):
        plt.figure(figsize=(12, 10))
        data_to_plot = self.interaction_grouped.head(self.top_interactions)
        sns.barplot(x='prevalence', y='interaction', data=data_to_plot, palette='viridis')
        plt.title(f'Top {len(data_to_plot)} Most Prevalent Interactions Between Chains')
        plt.xlabel('Average Prevalence (%)')
        plt.ylabel('Interaction (Residue-Residue)')
        plt.tight_layout()

        output_png = f'prev_interaction_top_{len(data_to_plot)}.png'
        plt.savefig(output_png, bbox_inches='tight')
        print(f"\nGráfico salvo como '{output_png}'")
        plt.show()

    def run(self):
        self.load_data()
        self.process_data()
        self.save_csv()
        self.plot_interactions()


def main():
    ProcessInteractions.main()
    parser = argparse.ArgumentParser(description='Analyze and visualize protein chain interactions.')
    parser.add_argument('--input', type=str, default='./detailed_interactions.tsv', help='Path to the input TSV file')
    parser.add_argument('--top', type=int, default=30, choices=[10, 20, 30, 40, 50], help='Number of top interactions to visualize')
    args = parser.parse_args()

    view_interactions = ViewInteractions(args.input, args.top)
    view_interactions.run()

#if __name__ == '__main__':
#    main()
