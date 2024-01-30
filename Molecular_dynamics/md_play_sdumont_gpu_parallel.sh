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
