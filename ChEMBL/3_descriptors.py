import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors

class MolecularDescriptors:
    def __init__(self, data_path):
        self.data_path = data_path
        self.descriptor_names = ['MW', 'LogP', 'HBD', 'HBA', 'TPSA', 'NRB']
        self.descriptor_data = {name: [] for name in self.descriptor_names}
        self.data = pd.read_csv(data_path, sep='\t')

    def calculate_descriptors(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return {
                    'MW': Descriptors.MolWt(mol),
                    'LogP': Descriptors.MolLogP(mol),
                    'HBD': Descriptors.NumHDonors(mol),
                    'HBA': Descriptors.NumHAcceptors(mol),
                    'TPSA': Descriptors.TPSA(mol),
                    'NRB': Descriptors.NumRotatableBonds(mol)
                }
        except Exception as e:
            print(f"Erro ao processar SMILES: {smiles}: {e}")
        return None

    def compute_descriptors(self):
        for smiles in self.data['salt_free_smiles']:
            descriptors = self.calculate_descriptors(smiles)
            if descriptors:
                for name in self.descriptor_names:
                    self.descriptor_data[name].append(descriptors[name])

    def save_descriptors(self, output_path):
        descriptor_df = pd.DataFrame(self.descriptor_data)
        descriptor_df.to_csv(output_path, sep='\t', index=False)

    def plot_histograms(self, additional_data_path=None, output_path=None):
        fig, axs = plt.subplots(3, 2, figsize=(13, 13)) # 3 linhas, 2 colunas
        axs = axs.flatten() # Transformar o array 2D de eixos em 1D

        # Criar eixos invisíveis para os rótulos "Density"
        density_axis = fig.add_subplot(111, frameon=False)  # Adicionar eixo que abrange toda a figura
        density_axis.set_xticks([])  # Desativar ticks do eixo x
        density_axis.set_yticks([])  # Desativar ticks do eixo y
        density_axis.grid(False)     # Sem grade
        density_axis.set_ylabel('Density', labelpad=40)  # Definir rótulo y

        for i, desc in enumerate(self.descriptor_names):
            axs[i].hist(self.descriptor_data[desc], bins=30, alpha=0.5, label='nr_ChEMBL', edgecolor='black', density=True)
            if additional_data_path:
                additional_data = pd.read_csv(additional_data_path, sep='\t')
                axs[i].hist(additional_data[desc], bins=30, alpha=0.5, label='PKIDB', edgecolor='black', density=True)
            axs[i].set_title(f'{desc}', fontsize=9)

        # Posicionar a legenda fora do plot no canto superior direito do último gráfico da primeira linha
        axs[1].legend(loc='upper right', bbox_to_anchor=(1.3, 1))

        plt.tight_layout()
        if output_path:
            plt.savefig(output_path)
        plt.show()
            #
            #plt.show()
       # else:
       #     self.plot_simple_histograms()
    
    def plot_simple_histograms(self):
        fig, axs = plt.subplots(3, 2, figsize=(13, 13)) # 3 linhas, 2 colunas
        axs = axs.flatten() # Transformar o array 2D de eixos em 1D para facilitar o acesso

        for i, name in enumerate(self.descriptor_names):
            axs[i].hist(self.descriptor_data[name], bins=30, edgecolor='black')
            axs[i].set_title(f'Histograma do Descritor: {name}', fontsize=9)
            axs[i].set_xlabel(name)
            axs[i].set_ylabel('Frequência')

        plt.tight_layout()
        plt.show()


def main():

    # Caminhos dos arquivos
    data_file_path = '../1_database/positive.tsv'
    additional_data_file_path = '../1_database/pkidb_2023-06-30.tsv'
    output_file_path = './chembl_nr_pkidb_descriptors.tsv'
    histogram_output_path = './nr_chembl_pkidb_descriptors.png'

    # Criar instância da classe
    molecular_descriptors = MolecularDescriptors(data_file_path)

    # Calcular descritores
    molecular_descriptors.compute_descriptors()

    # Salvar descritores calculados
    molecular_descriptors.save_descriptors(output_file_path)

    # Plotar histogramas normalizados e sobrepostos
    molecular_descriptors.plot_histograms(additional_data_file_path, histogram_output_path)

if __name__ == '__main__':
    main()
