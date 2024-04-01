#!/bin/bash
#SBATCH --nodes=1                  # Número de nós
#SBATCH -p sequana_dockvs          # Fila (partition) a ser utilizada
#SBATCH -J mdynamics               # Nome do job
#SBATCH --exclusive                # Utilização exclusiva dos nós

# Este script automatiza a execução de simulações de dinâmica molecular utilizando Amber em múltiplas réplicas.
# Ele distribui as tarefas de forma eficiente entre 4 GPUs disponíveis e é projetado para ser executado no supercomputador SDumont.
# Certifique-se de que o ambiente Amber está corretamente configurado e que os diretórios de réplicas estão preparados conforme esperado.

# INPUTS:

working_dir="/scratch/dockvs/leon.costa/way/run_11"

total_replicas=10                                        # Número total de réplicas
num_gpus=4                                               # Número de GPUs disponíveis

export PRMTOP="$working_dir/5cc8_way_run_11_tleap.prmtop"   # TOPOLOGY
export RST7="$working_dir/5cc8_way_run_11_tleap.rst7"       # COORDINATE + VELOCITY
export DO_CUDA="pmemd.cuda"                              # GPU
export DO_PARALLEL="mpirun -np 48 sander.MPI"            # CPU

source /scratch/dockvs/softwares/amber22/app/amber.sh    # amber path

run_replica() {
    local REPLICA=$1
    export CUDA_VISIBLE_DEVICES=$(( (REPLICA - 1) % $num_gpus ))

    if [ ! -d "$working_dir/$REPLICA" ]; then
        echo "Diretório $working_dir/$REPLICA não encontrado!"
        exit 1
    else
        cd $working_dir/$REPLICA/
    fi

    # Arrays para inputs, restarts, trajectories e outputs
    inputs=("0_minimization.in" "1_relax-part1.in" "2_relax-part2.in" "3_relax-part3.in" "4_relax-part4.in" "5_pre_prod.in" "6_production.in")
    rst_in=("$RST7" "minimization.rst" "relax-part1.rst" "relax-part2.rst" "relax-part3.rst" "relax-part4.rst" "pre-prod.rst")
    rst_out=("minimization.rst" "relax-part1.rst" "relax-part2.rst" "relax-part3.rst" "relax-part4.rst" "pre-prod.rst" "production.rst")
    trajectories=("minimization.crd" "relax-part1.crd" "relax-part2.crd" "relax-part3.crd" "relax-part4.crd" "pre-prod.crd" "production.crd")
    outputs=("minimization.out" "relax-part1.out" "relax-part2.out" "relax-part3.out" "relax-part4.out" "pre-prod.out" "production.out")

    for i in {0..6}; do
        input=${inputs[$i]}
        ref_rst=${rst_in[$i]}
        out_rst=${rst_out[$i]}
        crd=${trajectories[$i]}
        output=${outputs[$i]}

        echo "Réplica $REPLICA, Etapa $((i)): $input"
        if [ $i -eq 0 ]; then
            $DO_CUDA -O -i $input -p $PRMTOP -c $ref_rst -r $out_rst -x $crd -o $output
        elif [ $i -ge 1 ] && [ $i -le 4 ]; then
            $DO_CUDA -O -i $input -p $PRMTOP -c $ref_rst -ref $ref_rst -r $out_rst -x $crd -o $output
        elif [ $i -eq 5 ]; then
            $DO_CUDA -O -i $input -p $PRMTOP -c $ref_rst -ref $ref_rst -r $out_rst -x $crd -o $output
        else
            $DO_CUDA -O -i $input -p $PRMTOP -c $ref_rst -ref $ref_rst -r $out_rst -x $crd -o $output
        fi

        echo "Réplica $REPLICA, Etapa $((i)) finalizada: $output"
    done

    echo "Réplica $REPLICA finalizada."
}

main() {
    for (( i=1; i<=total_replicas; i+=num_gpus )); do
        for (( j=i; j<i+num_gpus && j<=total_replicas; j++ )); do
            run_replica $j &
        done
        wait
    done
    echo 'Processamento das réplicas finalizado.'
}
main
