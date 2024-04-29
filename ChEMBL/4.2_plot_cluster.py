import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns

# Carregar dados
data = pd.read_csv('./chembl_cluster_hits.tsv', sep='\t')

# Converter SMILES em fingerprints
def smiles_to_fp(smiles, radius=2, nBits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
    else:
        return [0]*nBits

# Aplicar a conversão e guardar os resultados em uma nova coluna
data['fingerprint'] = data['canonical_smiles'].apply(smiles_to_fp)

# Preparar os dados para o t-SNE
fp_list = list(data['fingerprint'])
fp_array = np.array([list(fp) for fp in fp_list])

# Executar t-SNE
tsne = TSNE(n_components=2, random_state=42)
tsne_results = tsne.fit_transform(fp_array)

# Incluir resultados t-SNE no DataFrame
data['tsne-2d-one'] = tsne_results[:,0]
data['tsne-2d-two'] = tsne_results[:,1]

# Definir uma paleta de cores personalizada
kinase_groups = data['kinase_group'].unique()
colors = plt.cm.get_cmap('tab20', len(kinase_groups))  # Usando o mapa de cores 'tab20' para diversidade

# Mapear cada grupo de quinase para uma cor específica
color_map = {group: colors(i) for i, group in enumerate(kinase_groups)}

# Plotar o gráfico t-SNE colorido pelos grupos de kinase
plt.figure(figsize=(16,10))
sns.scatterplot(
    x="tsne-2d-one", y="tsne-2d-two",
    hue="kinase_group",
    palette=color_map,  # Usando a paleta de cores personalizada
    data=data,
    legend="full",
    alpha=0.5
)
plt.title('t-SNE colored by Kinase Groups')
plt.show()
