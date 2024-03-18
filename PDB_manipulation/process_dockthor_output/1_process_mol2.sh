#!/bin/bash

# Caminho para o arquivo PDB da proteína
protein_pdb="../KP02043_5cc8_maestro_apo.pdb"

global_index=1

# Cria o diretório de saída se ele não existir
mkdir -p mol2_to_pdb

# Ordenar e processar os arquivos .mol2 garantindo a sequência correta para cada run
find . -maxdepth 1 -name '*.mol2' -print0 | sort -zV | while IFS= read -r -d '' mol2file; do
    base=$(basename "$mol2file" .mol2)

    # Dividir o nome do arquivo para extrair as partes relevantes
    IFS='_' read -ra ADDR <<< "$base"
    docked_run="${ADDR[4]}"
    molecule="${ADDR[6]}"

    # Gerar um nome base que reflete a ordem correta
    new_base="${ADDR[0]}_${ADDR[1]}_${ADDR[2]}_${ADDR[3]}_${docked_run}_molecule_${molecule}"

    # Usar o Open Babel para converter cada arquivo .mol2 para .pdb
    obabel "$mol2file" -O "mol2_to_pdb/${global_index}_${new_base}.pdb" -m

    ((global_index++))
done

# A parte de criação do complexo proteína-ligante permanece a mesma
mkdir -p pdb_complexed

for ligand_pdb in mol2_to_pdb/*.pdb; do
    base_name=$(basename "$ligand_pdb" .pdb)

    complex_pdb_name="${base_name}_complexed.pdb"
    complex_pdb_path="pdb_complexed/$complex_pdb_name"

    protein_content=$(cat "$protein_pdb")
    ligand_content=$(cat "$ligand_pdb")

    if ! echo "$protein_content" | grep -q "TER"; then
        protein_content+="TER\n"
    fi

    echo -e "$protein_content\n$ligand_content\nEND" > "$complex_pdb_path"
done
