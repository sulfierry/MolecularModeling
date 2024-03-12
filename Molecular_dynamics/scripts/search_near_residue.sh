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
