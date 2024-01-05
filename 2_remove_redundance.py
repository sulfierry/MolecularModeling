import pandas as pd
from rdkit import Chem
from rdkit.Chem import SaltRemover

def remove_salts_and_identify(smiles):
    remover = SaltRemover.SaltRemover()
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        salt_free_mol = remover.StripMol(mol)
        return Chem.MolToSmiles(salt_free_mol), salt_free_mol != mol
    return smiles, False
