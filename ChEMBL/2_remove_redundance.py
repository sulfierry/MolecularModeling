import pandas as pd
from rdkit import Chem
from rdkit.Chem import SaltRemover
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

class RemoveRedundance:

    def __init__(self, input_file_path, output_directory):
        self.input_file_path = input_file_path
        self.output_directory = output_directory
        self.remover = SaltRemover.SaltRemover()

    def remove_salts_and_identify(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            salt_free_mol = self.remover.StripMol(mol)
            salt_free_smiles = Chem.MolToSmiles(salt_free_mol)
            return salt_free_smiles, salt_free_smiles != smiles
        return smiles, False

    def process_smiles(self, smiles_list):
        with ThreadPoolExecutor(max_workers=4) as executor:
            results = list(tqdm(executor.map(self.remove_salts_and_identify, smiles_list), total=len(smiles_list), desc="Removing salts"))
        return results

    def process_data(self):
        data = pd.read_csv(self.input_file_path, sep='\t')
        results = self.process_smiles(data['canonical_smiles'].tolist())
        data['canonical_smiles'], data['contains_salt'] = zip(*results)
        return data

    def save_data(self, data, filename):
        data.to_csv(f"{self.output_directory}/{filename}", sep='\t', index=False)

    def execute(self):
        print("Processing data...")
        processed_data = self.process_data()

        print("Saving all data files...")
        self.save_data(processed_data, 'nr_kinase_all_compounds_with_salt.tsv')
        self.save_data(processed_data, 'nr_kinase_all_compounds_salt_free.tsv')
        self.save_data(processed_data[processed_data['standard_value'] <= 10000], 'positive.tsv')
        self.save_data(processed_data[processed_data['standard_value'] > 10000], 'negative.tsv')

        print("Data processing and file generation completed.")

if __name__ == "__main__":
    remover = RemoveRedundance('../0_database/kinase_all_compounds.tsv', '.')
    remover.execute()
