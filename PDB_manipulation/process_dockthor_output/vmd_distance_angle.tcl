#!/bin/bash

# Diretório onde estão os arquivos .pdb

PDB_DIR="4_pdb_frame_formated"
tcl_scipt="angle.tcl"
input_file="angle.dat"
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
awk '{print NR " " $2}' "$input_file" > "$output_file"
echo "Arquivo formatado: $output_file"
