#!/bin/sh
# No Santos Dummont
source /scratch/dockvs/softwares/amber22/app/amber.sh

# Configurações gerais
do_cuda="pmemd.cuda"
do_parallel="mpirun -np 48 sander.MPI"
prmtop="5cc8.prmtop"

# Arrays para inputs, restarts, trajectories e outputs
inputs=("0_minimization.in" "1_relax-part1.in" "2_relax-part2.in" "3_relax-part3.in" "4_relax-part4.in" "5_pre_prod.in" "6_production.in")
rst_in=("5cc8.rst7" "minimization.rst" "relax-part1.rst" "relax-part2.rst" "relax-part3.rst" "relax-part4.rst" "pre-prod.rst")
rst_out=("minimization.rst" "relax-part1.rst" "relax-part2.rst" "relax-part3.rst" "relax-part4.rst" "pre-prod.rst" "production.rst")
trajectories=("minimization.crd" "relax-part1.crd" "relax-part2.crd" "relax-part3.crd" "relax-part4.crd" "pre-prod.crd" "production.crd")
outputs=("minimization.out" "relax-part1.out" "relax-part2.out" "relax-part3.out" "relax-part4.out" "pre-prod.out" "production.out")


abort()
{
    echo >&2 '

************************
*** ABORTAR EXECUÇÃO ***
************************
'
    echo "Erro inesperado..." >&2
    exit 1
}

trap 'abort' 0

set -e

echo "Iniciando protocolo..."

