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

def redundance_remov(file_path):
    # Carrega os dados
    original_data = pd.read_csv(file_path, sep='\t')

    # Identifica se o SMILE possui algum sal
    original_data['contains_salt'] = original_data['canonical_smiles'].apply(lambda x: remove_salts_and_identify(x)[1])

    # Agrupa por 'canonical_smiles' e obtém lista de cinases-alvo
    smiles_kinases = original_data.groupby('canonical_smiles')['target_kinase'].apply(list).reset_index()
    smiles_kinases['num_of_kinases'] = smiles_kinases['target_kinase'].apply(lambda x: len(set(x)))

    # Faz merge com o DataFrame original para adicionar as colunas necessárias
    merged_data = pd.merge(original_data, smiles_kinases, on='canonical_smiles', how='left')

    # Renomeia a coluna 'target_kinase' para evitar conflitos após o merge
    merged_data.rename(columns={'target_kinase_x': 'target_kinase_list', 'target_kinase_y': 'target_kinase'}, inplace=True)

    # Seleciona a linha com o menor 'standard_value' para cada SMILE
    final_data = merged_data.sort_values(by=['canonical_smiles', 'standard_value']).drop_duplicates(subset='canonical_smiles', keep='first')

    # Adiciona as colunas 'kinase_total' e 'compound_total'
    final_data['kinase_total'] = None
    final_data['compound_total'] = None

    # Verifica se o DataFrame não está vazio antes de atribuir os valores
    if not final_data.empty:
        final_data.at[final_data.index[0], 'kinase_total'] = original_data['target_kinase'].nunique()
        final_data.at[final_data.index[0], 'compound_total'] = original_data['canonical_smiles'].nunique()

    # Reordena as colunas
    colunas = ['chembl_id', 'molregno', 'canonical_smiles', 'organism', 'target_kinase', 'num_of_kinases', 'contains_salt', 'standard_value', 'standard_type', 'pchembl_value', 'compound_name', 'kinase_total', 'compound_total']
    final_data = final_data[colunas]

    return final_data
    
