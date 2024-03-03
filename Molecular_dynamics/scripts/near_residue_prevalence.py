import re
from collections import defaultdict

# Caminho para o arquivo de dados
file_path = './residues_near_ligand.dat'

# Dicionário para contar a ocorrência de cada resíduo
residue_count = defaultdict(int)

# Processar o arquivo
with open(file_path, 'r') as file:
    for line in file:
        if line.startswith('    Residue '):
            # Extrair os resíduos da linha
            residues = re.findall(r'([A-Z]+ \d+)', line)
            for residue in residues:
                residue_count[residue] += 1

# Número total de frames
total_frames = max(int(line.split()[1].rstrip(':')) for line in open(file_path) if line.startswith("Frame")) + 1

# Calcular a prevalência de cada resíduo
residue_prevalence = {residue: count / total_frames for residue, count in residue_count.items()}

# Vamos converter os resultados em uma lista de tuplas para facilitar a visualização e escrita em um arquivo CSV
residue_prevalence_list = [(res.split()[1], res.split()[0], prevalence) for res, prevalence in residue_prevalence.items()]

# Ordenar a lista por prevalência em ordem decrescente
residue_prevalence_list.sort(key=lambda x: x[2], reverse=True)

# Imprimindo os resultados
for item in residue_prevalence_list:
    print(f"{item[0]:>10} {item[1]:>10} {item[2]*100:.2f}")

