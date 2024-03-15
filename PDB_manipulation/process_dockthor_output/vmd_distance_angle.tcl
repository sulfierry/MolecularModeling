vmd_angle_distance.sh#!/bin/bash


PDB_DIR="pdb_frames"
tcl_scipt="distance.tcl"
input_file="distance.dat"
output_file="distance_formatted.dat"

# Lista os arquivos .pdb na ordem correta usando o comando ls com a opção -v
for pdb_file in $(ls -v $PDB_DIR/*.pdb); do
    # Extrai o nome base do arquivo sem a extensão .pdb
    base_name=$(basename "$pdb_file" .pdb)
    echo "Processando $base_name..."

    # Executa o script VMD com o arquivo .pdb atual
    vmd -dispdev text -e $tcl_scipt "$pdb_file"
done

echo "Todos os arquivos foram processados."

# Usa awk para processar o arquivo e reescrever a primeira coluna com o número da linha
awk '{print NR " " $0}' "$input_file" > "$output_file"
echo "Arquivo formatado: $output_file"
