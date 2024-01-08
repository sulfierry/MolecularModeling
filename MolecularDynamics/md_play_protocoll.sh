#!/bin/sh
# source /scratch/dockvs/softwares/amber22/app/amber.sh

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

# Loop para executar as etapas
for i in {0..6}; do
    input=${inputs[$i]}
    ref_rst=${rst_in[$i]}
    out_rst=${rst_out[$i]}
    crd=${trajectories[$i]}
    output=${outputs[$i]}

    echo "Iniciando etapa $((i))..."

    if [ $i -eq 0 ]; then
        # Minimização
        echo "Iniciando minimização..."
        # 1000 passos de minimizacao energetica"
        $do_cuda -O -i $input -p $prmtop -c $ref_rst -r $out_rst -x $crd -o $output
        
    elif [ $i -ge 1 ] && [ $i -le 4 ]; then
        # Relaxamento partes 1-4 (1 ns)
        # 1: 300ps de pré-relaxamento NPT com restrição para proteína e ligante
        # 2: 300ps de relaxamento NPT com restrição para proteína
        # 3: 200ps de relaxamento NPT sem restrição para as cadeias laterais dos resíduos em até 5 Angstrons ao redor do ligante
        # 4: 200ps de relaxamento NPT sem restrição para os resíduos em até 5 Angstrons ao redor do ligante
        echo "Executando Relaxamento Parte $((i))"
        $do_cuda -O -i $input -p $prmtop -c $ref_rst -ref $ref_rst -r $out_rst -x $crd -o $output
        
    elif [ $i -eq 5 ]; then
        # Pré-produção (5 ns)
        echo "Executando Pré-produção..."
        $do_cuda -O -i $input -p $prmtop -c $ref_rst -ref $ref_rst -r $out_rst -x $crd -o $output
        
    else
        # Produção (100 ns)
        echo "Executando Produção..."
        $do_cuda -O -i $input -p $prmtop -c $ref_rst -ref $ref_rst -r $out_rst -x $crd -o $output
        
    fi

    echo "Etapa $((i)) finalizada com sucesso!"
done


trap : 0

echo >&2 '
*************
***  FIM  ***
*************
'
