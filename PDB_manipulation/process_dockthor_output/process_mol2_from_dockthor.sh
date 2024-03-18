#!/bin/bash

# Caminho para o arquivo PDB da proteína
protein_pdb="../KP02043_5cc8_maestro_apo.pdb"

# Cria o diretório de saída se ele não existir
mkdir -p mol2_to_pdb

# Ordenar e processar os arquivos .mol2 garantindo a sequência correta para cada run
find . -maxdepth 1 -name '*.mol2' -print0 | sort -zV | while IFS= read -r -d '' mol2file; do
    base=$(basename "$mol2file" .mol2)

    # Extrair corretamente o número da run
    if [[ "$base" =~ result-ligand_([0-9a-f]+)_([0-9]+)_docked_run_([0-9]+) ]]; then
        ligand_id="${BASH_REMATCH[1]}"
        attempt_number="${BASH_REMATCH[2]}"
        run_number="${BASH_REMATCH[3]}"

        # Gerar um nome base que reflete a ordem correta, incluindo o valor da "run"
        new_base="result-ligand_${ligand_id}_${attempt_number}_docked_run_${run_number}"

        # Usar o Open Babel para converter cada arquivo .mol2 para .pdb
        obabel "$mol2file" -O "mol2_to_pdb/${new_base}_molecule_.pdb" -m
    else
        echo "O nome do arquivo '$base' não corresponde ao padrão esperado."
    fi
done

# Pós-processamento para adicionar o índice corretamente de maneira crescente
global_index=1
for pdbfile in $(ls mol2_to_pdb/*.pdb | sort -V); do
    mv "$pdbfile" "mol2_to_pdb/${global_index}_$(basename "$pdbfile")"
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


mkdir -p pdb_frame_formated

# Percorre todos os arquivos .pdb na pasta "pdb_complexed", ordenados pelo número antes do primeiro underscore
for complexed_pdb in $(find pdb_complexed -name '*.pdb' | sed 's/^.*\/\([0-9]*\)_.*\.pdb$/\1 &/' | sort -n | cut -d' ' -f2-); do
    # Extrai o nome base do arquivo (sem a extensão .pdb)
    base_name=$(basename "$complexed_pdb" .pdb)

    # Define o caminho do arquivo de saída na pasta "pdb_frame_formated"
    formatted_pdb_path="pdb_frame_formated/${base_name}_formated.pdb"

    # Usa o pdb4amber para formatar o arquivo PDB e salva o resultado no arquivo de saída
    pdb4amber -i "$complexed_pdb" -o "$formatted_pdb_path"
done

mkdir -p pdb4amber_infos
mv ./pdb_frame_formated/*_nonprot.* pdb4amber_infos
mv ./pdb_frame_formated/*_formated_renum.* pdb4amber_infos
mv ./pdb_frame_formated/*formated_sslink* pdb4amber_infos
