import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import SaltRemover

class RemoveRedundance:
    def __init__(self, input_file_path, output_directory):
        self.input_file_path = input_file_path
        self.output_directory = output_directory
        self.remover = SaltRemover.SaltRemover()

    def remove_salts_and_identify(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            original_smiles = Chem.MolToSmiles(mol)
            salt_free_mol = self.remover.StripMol(mol)
            salt_free_smiles = Chem.MolToSmiles(salt_free_mol)
            return salt_free_smiles, original_smiles != salt_free_smiles
        return smiles, False

    def process_data(self):
        data = pd.read_csv(self.input_file_path, sep='\t')
        tqdm.pandas(desc="Removing salts")
        results = data['canonical_smiles'].progress_apply(self.remove_salts_and_identify)
        data['canonical_smiles'], data['contains_salt'] = zip(*results)

        return data

    def create_salt_free_output(self, data):
        tqdm.pandas(desc="Processing to salt-free SMILES")
        data['canonical_smiles'] = data['canonical_smiles'].progress_apply(lambda x: self.remove_salts_and_identify(x)[0])
        return data

    def save_data(self, data, filename):
        data.to_csv(f"{self.output_directory}/{filename}", sep='\t', index=False)

    def execute(self):
        print("Processing data...")
        processed_data = self.process_data()

        print("Creating salt-free output...")
        salt_free_data = self.create_salt_free_output(processed_data)

        print("Saving all data files...")
        self.save_data(processed_data, 'nr_kinase_all_compounds_with_salt.tsv')
        self.save_data(salt_free_data, 'nr_kinase_all_compounds_salt_free.tsv')
        self.save_data(salt_free_data[salt_free_data['standard_value'] <= 10000], 'positive.tsv')
        self.save_data(salt_free_data[salt_free_data['standard_value'] > 10000], 'negative.tsv')

        print("Data processing and file generation completed.")


class FormatChembl:
    def __init__(self, input_file_path, pkidb_file, output_file_path):
        self.input_file_path = input_file_path
        self.pkidb_file = pkidb_file
        self.output_file_path = output_file_path

    @staticmethod
    def process_kinase_group(group_name):
        special_cases = {
            "Enzyme": "Other", "Protein": "Other", "Other": "Other", "Unclassified": "Other",
            "Atypical": "Other", "Kinase": "Other", "Nuclear e Other nuclear protein": "Nuclear protein",
            "Membrane": "Membrane receptor", "Voltage-gated potassium channel": "Potassium channel",
            "Tudor domain": "Tudor domain", "Bromodomain": "Bromodomain",
            "Transcription factor": "Transcription factor", "Plant homeodomain": "Plant homeodomain",
            "Protein Kinase": "Protein Kinase"
        }
        for key, value in special_cases.items():
            if group_name.startswith(key):
                return value
        return group_name.split()[0] if group_name else group_name

    def filter_pkidb_chembl_ids(self, input_data):
        # Load PKIDB data
        pkidb_data = pd.read_csv(self.pkidb_file, sep='\t')
        # Extract CHEMBL_IDs from PKIDB
        pkidb_ids = set(pkidb_data['CHEMBL_ID'])
        # Filter out entries in input_data that are in PKIDB
        return input_data[~input_data['chembl_id'].isin(pkidb_ids)]

    def format_data(self):
        # Load input data
        data = pd.read_csv(self.input_file_path, sep='\t')
        # Process kinase group names
        if 'kinase_group' in data.columns:
            data['kinase_group'] = data['kinase_group'].apply(self.process_kinase_group)
        # Filter out CHEMBL_IDs from PKIDB
        data = self.filter_pkidb_chembl_ids(data)
        # Save formatted data
        data.to_csv(self.output_file_path, sep='\t', index=False)
        print("Processed data saved successfully.")

def format_chembl_file():
    file_path = '../0_database/kinase_all_compounds.tsv'
    pkidb_file = '../0_database/pkidb/pkidb_2024-03-18.tsv'
    output_file_path = '../0_database/kinase_all_compounds_formatted.tsv'
    formatter = FormatChembl(file_path, pkidb_file, output_file_path)
    formatter.format_data()

if __name__ == "__main__":
    format_chembl_file()
    remover = RemoveRedundance('../0_database/kinase_all_compounds_formatted.tsv', '.')
    remover.execute()
