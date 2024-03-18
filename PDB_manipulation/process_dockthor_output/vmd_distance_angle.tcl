#!/bin/bash

# Diretório onde estão os arquivos .pdb
PDB_DIR="4_pdb_frame_formated"
tcl_script="angle.tcl"
input_file="angle_ver2.dat"
output_file="angle_formatted_ver2.dat"

export PDB_DIR
export tcl_script

# Define a função a ser executada pelo parallel
process_pdb() {
    pdb_file="$1"
    base_name=$(basename "$pdb_file" .pdb)
    echo "Processando $base_name..."

    # Executa o script VMD com o arquivo .pdb atual
    vmd -dispdev text -e $tcl_script "$pdb_file"
}
export -f process_pdb

# Lista os arquivos .pdb na ordem correta e os processa em paralelo
find $PDB_DIR -name '*.pdb' | sort -V | parallel process_pdb

echo "Todos os arquivos foram processados."

# Usa awk para processar o arquivo e reescrever a primeira coluna com o número da linha
awk '{print NR " " $2}' "$input_file" > "$output_file"
echo "Arquivo formatado: $output_file"
