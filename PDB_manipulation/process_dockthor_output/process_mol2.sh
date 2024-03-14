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


