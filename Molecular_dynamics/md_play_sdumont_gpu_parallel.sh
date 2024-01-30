#!/bin/bash
#SBATCH --nodes=1                  # Número de nós
#SBATCH -p sequana_dockvs          # Fila (partition) a ser utilizada
#SBATCH -J molecular_dynamics      # Nome do job
#SBATCH --exclusive                # Utilização exclusiva dos nós


# Carregar o ambiente necessário
source /scratch/dockvs/softwares/amber22/app/amber.sh

# Configurações gerais
export DO_CUDA="pmemd.cuda"
export DO_PARALLEL="mpirun -np 48 sander.MPI"
export PRMTOP="5cc8.prmtop"
export RST7="5cc8.rst7"


# Função para executar uma réplica
run_replica() {
    local REPLICA=$1

    # Define a GPU a ser utilizada
    export CUDA_VISIBLE_DEVICES=$(( (REPLICA - 1) % 4 ))

    # cd /scratch/dockvs/leon.costa/md_thil_10replicates_100ns/$REPLICA/
    cd /scratch/dockvs/leon.costa/md_thil_10replicates_100ns/2_replica/$REPLICA/

    # Arrays para inputs, restarts, trajectories e outputs
    inputs=("0_minimization.in" "1_relax-part1.in" "2_relax-part2.in" "3_relax-part3.in" "4_relax-part4.in" "5_pre_prod.in" "6_production.in")
    rst_in=("$RST7" "minimization.rst" "relax-part1.rst" "relax-part2.rst" "relax-part3.rst" "relax-part4.rst" "pre-prod.rst")
    rst_out=("minimization.rst" "relax-part1.rst" "relax-part2.rst" "relax-part3.rst" "relax-part4.rst" "pre-prod.rst" "production.rst")
    trajectories=("minimization.crd" "relax-part1.crd" "relax-part2.crd" "relax-part3.crd" "relax-part4.crd" "pre-prod.crd" "production.crd")
    outputs=("minimization.out" "relax-part1.out" "relax-part2.out" "relax-part3.out" "relax-part4.out" "pre-prod.out" "production.out")

    # Loop para executar as etapas de MD para a réplica especificada
    for i in {0..6}; do
        input=${inputs[$i]}
        ref_rst=${rst_in[$i]}
        out_rst=${rst_out[$i]}
        crd=${trajectories[$i]}
        output=${outputs[$i]}

        echo "Réplica $REPLICA, Etapa $((i)): $input"

        if [ $i -eq 0 ]; then
            # Minimização
            $DO_CUDA -O -i $input -p $PRMTOP -c $ref_rst -r $out_rst -x $crd -o $output
        elif [ $i -ge 1 ] && [ $i -le 4 ]; then
            # Relaxamento partes 1-4
            $DO_CUDA -O -i $input -p $PRMTOP -c $ref_rst -ref $ref_rst -r $out_rst -x $crd -o $output
        elif [ $i -eq 5 ]; then
            # Pré-produção
            $DO_CUDA -O -i $input -p $PRMTOP -c $ref_rst -ref $ref_rst -r $out_rst -x $crd -o $output
        else
            # Produção
            $DO_CUDA -O -i $input -p $PRMTOP -c $ref_rst -ref $ref_rst -r $out_rst -x $crd -o $output
        fi

        echo "Réplica $REPLICA, Etapa $((i)) finalizada: $output"
    done

    echo "Réplica $REPLICA finalizada."
}
