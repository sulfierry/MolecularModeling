#!/bin/bash

main() {
    # Define o caminho absoluto onde as pastas estão localizadas e outras variáveis
    export ROOT_DIR="/media/leon/FEDF-FDB3/way_analogues/01_gmmsb/run_15/water_remov"
    export base_name="gmmsb01_run15_wr" # os arquivos precisar ter _{i} apos o $base_name
    export contact_threshold_percentage="10"
    export ligand_distance="3.00"
    export ligand_name="LIG"
    export start_folder="1"
    export end_folder="10"

    # Chama a função que executa o protocolo
    execute_protocol
}

execute_protocol() {

    # Gera o script TCL antes de iniciar o loop
    generate_tcl_script
    
    echo "Iniciando protocolo..."

    # Loop pelas pastas
    for ((i=start_folder; i<=end_folder; i++)); do
        echo "Processando pasta $i..."

        # Verifica se a pasta existe antes de tentar acessá-la
        if [ -d "$ROOT_DIR/$i" ]; then
            # Acessa a pasta
            cd "$ROOT_DIR/$i"
            
            echo "Diretório atual: $(pwd)"
            # ls -l

            # Define os nomes dos arquivos de topologia e trajetória baseados na pasta atual
            TOPOLOGY_FILENAME="$base_name"_"${i}.prmtop"
            TRAJECTORY_FILENAME="$base_name"_"${i}.dcd"

            # Executa o VMD com o script TCL, passando os nomes dos arquivos como argumentos, se necessário
            vmd -dispdev text -e ../search_near_residue.tcl $TOPOLOGY_FILENAME $TRAJECTORY_FILENAME
            
            # Retorna ao diretório raiz antes de continuar o loop
            cd "$ROOT_DIR"
        else
            echo "A pasta $ROOT_DIR/$i não existe."
        fi
    done

    echo "Processamento completo."
    # Remove o script TCL após o uso
    rm search_near_residue.tcl
    
    near_residue_prevalence
    # Executa o script Python para análise de prevalência dos resíduos
    if [ -f "near_residue_prevalence.py" ]; then
        python near_residue_prevalence.py
        # Remove o script Python após o uso, se necessário
        rm near_residue_prevalence.py
    else
        echo "O arquivo near_residue_prevalence.py não foi encontrado."
    fi
}




# Função para gerar o script search_near_residue.tcl
generate_tcl_script() {
    cat << EOF > search_near_residue.tcl
# vmd -dispdev text -e script.tcl 
# Find the number of atoms within 2-3 A of the ligand. 

# Define o nome do ligante
set ligand_name $ligand_name

# Cria uma seleção de todos os resíduos de proteína próximos ao ligante
set closer_residue [atomselect 0 "protein within $ligand_distance of resname \$ligand_name"]

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


