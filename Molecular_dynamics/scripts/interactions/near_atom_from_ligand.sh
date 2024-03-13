#!/bin/bash

# Configuração inicial
ROOT_DIR="$(pwd)"
BASE_NAME="gmmsb01_run15_wr_10"
TOPOLOGY="$BASE_NAME.prmtop"
COORDINATES="$BASE_NAME.dcd"
TOTAL_FRAMES=2000
SLICE=1000
LIGAND_NAME="LIG"
CUTOFF_DISTANCE="3.0"

# Função para gerar o script TCL dinamicamente
generate_tcl_script() {
    local start_frame=$1
    local end_frame=$2
    cat << EOF
# Verifica e cria a pasta "vmd_frames" se ela não existir
if {![file exists vmd_frames]} {
    file mkdir vmd_frames
}

# Defina o nome do ligante
set ligand_name "${LIGAND_NAME}"

# Defina a distância de corte
set cutoff_distance ${CUTOFF_DISTANCE}

# Obtém o número total de quadros na trajetória
set num_frames [molinfo 0 get numframes]

# Itera sobre cada quadro da simulação dentro do bloco definido
for {set frame_number $start_frame} {\$frame_number <= $end_frame} {incr frame_number} {
    # Atualiza para o quadro atual
    molinfo 0 set frame \$frame_number

    # Abre o arquivo para escrita para o quadro atual dentro da pasta "vmd_frames"
    set outfile [open "vmd_frames/detailed_interactions_frame_\$frame_number.dat" w]

    # Cria uma seleção de todos os átomos do ligante
    set ligand_atoms [atomselect 0 "resname \$ligand_name"]

    # Obtém informações detalhadas dos átomos do ligante
    foreach ligand_atom [\$ligand_atoms get {index name}] {
        set ligand_index [lindex \$ligand_atom 0]
        set ligand_atomname [lindex \$ligand_atom 1]
        
        # Cria uma seleção de todos os átomos de cadeias laterais de aminoácidos próximos a este átomo do ligante
        set sidechain_atoms [atomselect 0 "sidechain within \$cutoff_distance of index \$ligand_index"]
        set sidechain_atoms_info [\$sidechain_atoms get {index element resname resid name}]
        
        foreach sidechain_atom \$sidechain_atoms_info {
            set sidechain_index [lindex \$sidechain_atom 0]
            set sidechain_element [lindex \$sidechain_atom 1]
            set sidechain_resname [lindex \$sidechain_atom 2]
            set sidechain_resid [lindex \$sidechain_atom 3]
            set sidechain_atomname [lindex \$sidechain_atom 4]
            
            # Calcula a distância entre o átomo do ligante e o átomo do aminoácido
            set distance [measure bond [list \$ligand_index \$sidechain_index]]

            # Escreve as informações no arquivo, incluindo a distância
            puts \$outfile "Frame \$frame_number: Ligand atom \$ligand_atomname index \$ligand_index interacts with \$sidechain_resname \$sidechain_resid \$sidechain_atomname Distance: \$distance"
        }
    }

    # Fecha o arquivo
    close \$outfile
}

# Encerra o script VMD
quit
EOF
}

# Processa os frames em blocos
for ((start=0; start<TOTAL_FRAMES; start+=SLICE)); do
    end=$((start + SLICE - 1))
    if [ $end -ge $TOTAL_FRAMES ]; then
        end=$((TOTAL_FRAMES - 1))
    fi

    # Gera o script TCL para o bloco atual de frames
    tcl_script=$(generate_tcl_script $start $end)

    # Salva o script TCL gerado em um arquivo temporário
    tcl_script_file="${ROOT_DIR}/temp_tcl_script.tcl"
    echo "$tcl_script" > "$tcl_script_file"

    # Executa o VMD com o script TCL gerado
    vmd -dispdev text -e "$tcl_script_file" "$TOPOLOGY" "$COORDINATES"

    # Limpa o arquivo temporário do script TCL
    rm "$tcl_script_file"
done
