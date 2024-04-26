import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import os
import seaborn as sns


class MolecularDescriptors:

    def __init__(self, data_path, batch_size=1024):
        self.data_path = data_path
        self.batch_size = batch_size
        self.descriptor_names = ['MW', 'LogP', 'HBD', 'HBA', 'TPSA', 'NRB']
        self.data = pd.read_csv(data_path, sep='\t')

    def calculate_descriptors(self, smiles_list):
        results = []
        for smiles in smiles_list:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    results.append({
                        'canonical_smiles': smiles,
                        'MW': Descriptors.MolWt(mol),
                        'LogP': Descriptors.MolLogP(mol),
                        'HBD': Descriptors.NumHDonors(mol),
                        'HBA': Descriptors.NumHAcceptors(mol),
                        'TPSA': Descriptors.TPSA(mol),
                        'NRB': Descriptors.NumRotatableBonds(mol)
                    })
            except Exception as e:
                print(f"Erro ao processar SMILES: {smiles}: {e}")
        return results

    def compute_descriptors(self):
        smiles_data = self.data['canonical_smiles'].dropna().unique()
        chunks = [smiles_data[i:i + self.batch_size] for i in range(0, len(smiles_data), self.batch_size)]
        results = []

        with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
            for chunk_result in tqdm(executor.map(self.calculate_descriptors, chunks), total=len(chunks), desc="Calculating Descriptors"):
                results.extend(chunk_result)

        valid_results = [res for res in results if res is not None]
        descriptors_df = pd.DataFrame(valid_results)
        self.descriptor_data = pd.merge(descriptors_df, self.data[['canonical_smiles', 'kinase_group', 'count_kinase_group']].drop_duplicates(), on='canonical_smiles', how='left')

    def save_descriptors(self, output_path):
        self.descriptor_data.to_csv(output_path, sep='\t', index=False)

    def plot_histograms(self, additional_data_path=None, output_path=None):
        fig, axs = plt.subplots(3, 2, figsize=(13, 13))
        axs = axs.flatten()
        density_axis = fig.add_subplot(111, frameon=False)
        density_axis.set_xticks([])
        density_axis.set_yticks([])
        density_axis.grid(False)
        density_axis.set_ylabel('Density', labelpad=40)

        for i, desc in enumerate(self.descriptor_names):
            axs[i].hist(self.descriptor_data[desc], bins=30, alpha=0.5, label='nr_ChEMBL', edgecolor='black', density=True)
            if additional_data_path:
                additional_data = pd.read_csv(additional_data_path, sep='\t')
                axs[i].hist(additional_data[desc], bins=30, alpha=0.5, label='PKIDB', edgecolor='black', density=True)
            axs[i].set_title(f'{desc}', fontsize=9)
        axs[1].legend(loc='upper right', bbox_to_anchor=(1.3, 1))
        plt.tight_layout()
        if output_path:
            plt.savefig(output_path)
        plt.show()


    def violin_plot(self):
        columns = ['MW', 'LogP', 'HBD', 'HBA', 'TPSA', 'NRB']
        fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(15, 10))
        colors = sns.color_palette("muted", len(columns))

        for ax, col, color in zip(axs.flat, columns, colors):
            q_low = self.data[col].quantile(0.01)
            q_high = self.data[col].quantile(0.99)
            range_expand = (q_high - q_low) * 0.10
            if col in ['MW', 'LogP', 'HBA', 'TPSA']:
                range_expand *= 1.5
            filtered_data = self.data[(self.data[col] >= q_low - range_expand) & (self.data[col] <= q_high + range_expand)][col]

            sns.violinplot(data=filtered_data, ax=ax, color=color, inner="quartile", cut=0)
            ax.set_ylabel(col)
            ax.set_xlabel('')

        plt.tight_layout()
        plt.savefig('./violin_plot.png')
        plt.show()


def main():
    data_file_path = '../1_remove_redundance/nr_kinase_all_compounds_salt_free.tsv'
    additional_data_file_path = '../0_database/pkidb/pkidb_2024-03-18.tsv'
    output_file_path = './chembl_nr_pkidb_descriptors.tsv'
    histogram_output_path = './nr_chembl_pkidb_descriptors.png'

    #molecular_descriptors = MolecularDescriptors(data_file_path, batch_size=1024)  # Ajuste o tamanho do batch conforme necessário
    #molecular_descriptors.compute_descriptors()
    #molecular_descriptors.save_descriptors(output_file_path)
    #molecular_descriptors.plot_histograms(additional_data_file_path, histogram_output_path)

    descriptors = MolecularDescriptors(output_file_path)
    descriptors.violin_plot()

if __name__ == '__main__':
    main()
