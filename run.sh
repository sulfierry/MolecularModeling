#!/bin/bash
#SBATCH --nodes=1                   # Numero de nós 21 fila nvidia)
#SBATCH -p sequana_dockvs           # Fila (partition) a ser utilizada
#SBATCH -J molecular_dynamics       # Nome job
#SBATCH --exclusive                 # Utilização exclusiva dos nós durante a execução do job

# squeue -a -u $USER   - Ver jobs do usuario
# squeue -a -p sequana_dockvs  - Ver jobs da fila
# scancel jobid
# salloc -p sequana_dockvs -J job_name --exclusive
# ssh node_returned
# sbatch run.sh

# exibe os nos alocados para o Job
echo $SLURM_JOB_NODELIST
nodeset -e $SLURM_JOB_NODELIST

cd $SLURM_SUBMIT_DIR  # the directory from which sbatch was invoked

source /scratch/dockvs/softwares/amber22/app/amber.sh

# executa o script
for i in {1..10}; do
    cd /scratch/dockvs/leon.costa/md_thil_10replicates_100ns/$i/
    ./md_play.sh
done

