#!/bin/bash

# Define o caminho absoluto onde as pastas estão localizadas e outras variáveis
export ROOT_DIR="/media/leon/FEDF-FDB3/way_analogues/01_gmmsb/run_17/water_remov"
export BASE_NAME="gmmsb01_run17_wr" # os arquivos precisam ter _{i} após o $BASE_NAME
export LIGAND_DISTANCE="3.00"
export LIGAND_NAME="LIG"
export START_FOLDER="1"
export END_FOLDER="10"

execute_protocol() {
    echo "Iniciando protocolo..."

    for ((i=START_FOLDER; i<=END_FOLDER; i++)); do
        echo "Processando pasta $i..."

        if [ -d "$ROOT_DIR/$i" ]; then
            cd "$ROOT_DIR/$i"
            echo "Diretório atual: $(pwd)"

            TOPOLOGY_FILENAME="${BASE_NAME}_${i}.prmtop"
            TRAJECTORY_FILENAME="${BASE_NAME}_${i}.dcd"

            generate_tcl_script
            vmd -dispdev text -e ../search_near_residue.tcl -args $TOPOLOGY_FILENAME $TRAJECTORY_FILENAME
            rm ../search_near_residue.tcl

            cd "$ROOT_DIR"
        else
            echo "A pasta $ROOT_DIR/$i não existe."
        fi
    done

    echo "Processamento completo."
}

generate_tcl_script() {
    cat <<EOF >../search_near_residue.tcl
# Define o nome do ligante
set ligand_name "$LIGAND_NAME"

# Cria uma seleção de todos os resíduos de proteína próximos ao ligante
set closer_residue [atomselect 0 "protein within $LIGAND_DISTANCE of resname \$ligand_name"]

# Abre o arquivo para escrita
set outfile [open "residues_near_ligand.dat" w]

# Obtém o número total de quadros na trajetória
set num_frames [molinfo 0 get numframes]

# Itera sobre cada quadro da simulação
for {set frame_number 0} {\$frame_number < \$num_frames} {incr frame_number} {
    # Atualiza a seleção para o quadro atual
    \$closer_residue frame \$frame_number
    \$closer_residue update

    # Obtém os resíduos únicos próximos ao ligante neste quadro
    set residues [\$closer_residue get {resname resid}]
    set unique_residues [lsort -unique \$residues]

    # Escreve o número do quadro e o número de contatos únicos
    puts \$outfile "Frame \$frame_number: [llength \$unique_residues] unique residues near \$ligand_name"

    # Para cada resíduo único, imprime sua identificação
    foreach {resname resid} \$unique_residues {
        puts \$outfile "    Residue \$resname \$resid"
    }
}

# Fecha o arquivo
close \$outfile

# Encerra o script VMD
quit
EOF
}

execute_protocol
