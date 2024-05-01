import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
import concurrent.futures
import os

def load_data(filepath):
    data = pd.read_csv(filepath, sep='\t')
    return data

def smiles_to_fingerprints(smiles, radius=2, n_bits=2048):
    if pd.isna(smiles):
        return [0] * n_bits  # Returns a null fingerprint if SMILES is NaN
    mol = Chem.MolFromSmiles(str(smiles))
    if mol:
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius, n_bits).ToList()
    return [0] * n_bits

def convert_smiles_batch(smiles_list, radius=2, n_bits=2048):
    """ Convert a batch of SMILES to fingerprints using multiple processes. """
    with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        results = list(tqdm(executor.map(smiles_to_fingerprints, smiles_list,
                                         [radius] * len(smiles_list), [n_bits] * len(smiles_list)),
                           total=len(smiles_list), desc="Converting SMILES"))
    return results

def preprocess_data(data):
    batch_size = 524288
    fingerprints = []

    for i in tqdm(range(0, len(data), batch_size), desc="Processing batches"):
        batch = data.iloc[i:i + batch_size]
        batch_fingerprints = convert_smiles_batch(batch['canonical_smiles'].tolist())
        fingerprints.extend(batch_fingerprints)

    data['fingerprint'] = fingerprints

    return data

def split_data(data):
    # Convert list of fingerprint lists to numpy array for ML models
    X = np.array(data['fingerprint'].tolist())
    y = data['kinase_group'].values  # Use kinase_group as target

    # Splitting the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    return X_train, X_test, y_train, y_test

def main():
    # Load data
    filepath = 'test_200.tsv'  # Adjust as needed
    data = load_data(filepath)

    # Preprocess data
    data = preprocess_data(data)

    # Split data
    X_train, X_test, y_train, y_test = split_data(data)

    # Optionally save the split data
    np.savez('split_data.npz', X_train=X_train, X_test=X_test, y_train=y_train, y_test=y_test)

if __name__ == '__main__':
    main()
