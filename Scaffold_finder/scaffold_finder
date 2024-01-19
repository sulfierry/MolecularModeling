import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
import concurrent.futures
import os
from tqdm import tqdm

# Função para encontrar o scaffold de um SMILES
def get_scaffold(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    return Chem.MolToSmiles(scaffold) if scaffold else None

# Função para calcular a similaridade de Tanimoto
def calculate_similarity(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    fp1 = AllChem.GetMorganFingerprint(mol1, 2)
    fp2 = AllChem.GetMorganFingerprint(mol2, 2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

# Função para processar cada linha do DataFrame
def process_row(row):
    scaffold = get_scaffold(row['canonical_smiles'])
    similarity = calculate_similarity(row['canonical_smiles'], target_smile)
    return scaffold, similarity

# Carregar o arquivo
file_path = './kinase_all_compounds.tsv'
data = pd.read_csv(file_path, sep='\t')

# Smile alvo
target_smile = "NC(=O)CC(C(=O)O)N" # WAY_213613 sub-estrutura

# Processamento paralelo
cpu_count = os.cpu_count()
with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_count) as executor:
    results = list(tqdm(executor.map(process_row, data.to_dict('records')), total=len(data)))

# Adicionando os resultados ao DataFrame e reorganizando as colunas
data[['scaffold_smile', 'similarity_smiles']] = pd.DataFrame(results, index=data.index)
data = data[['chembl_id','molregno', 'target_kinase', 'canonical_smiles', 'scaffold_smile', 'similarity_smiles', 
             'standard_value', 'standard_type', 'pchembl_value', 'compound_name']]

# Filtrando por similaridade (exemplo: >= 10%)
data = data[data['similarity_smiles'] >= 0.10]

# Removendo a coluna 'kinase_total'
data.drop(columns=['kinase_total'], inplace=True, errors='ignore')

# Salvando o resultado em um novo arquivo .tsv
output_file_path = 'output_scaffolds_similarity.tsv'
data.to_csv(output_file_path, sep='\t', index=False)
