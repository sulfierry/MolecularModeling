#!/bin/bash

# Caminho para o arquivo PDB da proteína
protein_pdb="./KP02043_5cc8_maestro_apo.pdb"

global_index=1

# Cria o diretório de saída se ele não existir
mkdir -p mol2_to_pdb

# Encontra todos os arquivos .mol2 na pasta atual e os processa um por um
for mol2file in *.mol2; do
  # Extrai o nome base do arquivo (sem a extensão)
  base=$(basename "$mol2file" .mol2)

  # Usa o Open Babel para converter cada arquivo .mol2 para .pdb
  obabel "$mol2file" -O "mol2_to_pdb/${base}_molecule_.pdb" -m

  # Após a conversão, renomeia os arquivos gerados para incluir o índice global
  for pdbfile in mol2_to_pdb/${base}_molecule_*.pdb; do
    mv "$pdbfile" "mol2_to_pdb/${global_index}_${pdbfile#mol2_to_pdb/}"
    let global_index++
  done
done


# Verifica se a pasta "pdb_complexed" existe; se não, cria-a
mkdir -p pdb_complexed

# Percorre todos os arquivos .pdb na pasta "mol2_to_pdb"
for ligand_pdb in mol2_to_pdb/*.pdb; do
  # Extrai o nome base do arquivo (sem a extensão .pdb)
  base_name=$(basename "$ligand_pdb" .pdb)

  # Cria um novo nome de arquivo para o complexo proteína-ligante, adicionando "_complexed" antes da extensão .pdb
  complex_pdb_name="${base_name}_complexed.pdb"

  # Caminho completo do novo arquivo complexo na pasta "pdb_complexed"
  complex_pdb_path="pdb_complexed/$complex_pdb_name"

  # Lê o conteúdo do arquivo PDB da proteína
  protein_content=$(cat "$protein_pdb")

  # Verifica se o conteúdo da proteína já termina com um registro TER; se não, adiciona
  if ! echo "$protein_content" | grep -q "TER"; then
    protein_content+="TER\n"
  fi

  # Lê o conteúdo do arquivo PDB do ligante
  ligand_content=$(cat "$ligand_pdb")

  # Combina os conteúdos, adicionando o ligante ao final da proteína
  # e assegurando que o arquivo termina com "END"
  echo -e "$protein_content\n$ligand_content\nEND" > "$complex_pdb_path"
done


# Verifica se a pasta "pdb_complexed_formated" existe; se não, cria-a
mkdir -p pdb_complexed_formated


# Percorre todos os arquivos .pdb na pasta "pdb_complexed"
for complexed_pdb in pdb_complexed/*.pdb; do
    # Extrai o nome base do arquivo (sem a extensão .pdb)
    base_name=$(basename "$complexed_pdb" .pdb)

    # Define o caminho do arquivo de saída na pasta "pdb_complexed_formated"
    formatted_pdb_path="pdb_complexed_formated/${base_name}_formated.pdb"

    # Usa o pdb4amber para formatar o arquivo PDB e salva o resultado no arquivo de saída
    pdb4amber -i "$complexed_pdb" -o "$formatted_pdb_path"
done

