#!/bin/sh
source /scratch/dockvs/softwares/amber22/app/amber.sh


# MD_setup ###########################################################################################################################

# process
do_cuda="pmemd.cuda"
do_parallel="mpirun -np 48 sander.MPI"

# topology file
prmtop="5cc8.prmtop"

# coordinate file
inpcrd="5cc8.rst7"

# inputs
input_minimization="0_minimization.in"
input_relax_part_1="1_relax-part1.in"
input_relax_part_2="2_relax-part2.in"
input_relax_part_3="3_relax-part3.in"
input_relax_part_4="4_relax-part4.in"
input_pre_prod="5_pre_prod.in"
input_production="6_production.in"

# restart files
min_rst="minimization.rst"
relax_rst_1="relax-part1.rst"
relax_rst_2="relax-part2.rst"
relax_rst_3="relax-part3.rst"
relax_rst_4="relax-part4.rst"
pre_prod_rst="pre-prod.rst"
production_rst="production.rst"

# trajectory files
min_crd="minimization.crd"
relax_crd_1="relax-part1.crd"
relax_crd_2="relax-part2.crd"
relax_crd_3="relax-part3.crd"
relax_crd_4="relax-part4.crd"
pre_prod_crd="pre-prod.crd"
production_crd="production.crd"

# outputs
min_out="minimization.out"
relax_out_1="relax-part1.out"
relax_out_2="relax-part2.out"
relax_out_3="relax-part3.out"
relax_out_4="relax-part4.out"
pre_prod_out="pre-prod.out"
production_out="production.out"

#########################################################################################################################################

abort()
{
    echo >&2 '

************************
*** ABORTAR EXECUÇÂO ***
************************
'
    echo "Erro inesperado..." >&2
    exit 1
}

trap 'abort' 0

set -e

# MD_PROTOCOL ########################################################################################################################

# 1   ns relaxamento
# 5   ns pre-producao
# 100 ns producao

echo "Iniciando protocolo..."

#echo "Iniciando minimização:"
#echo "1000 passos de minimizacao energetica"
#$do_cuda -O -i $input_minimization -p $prmtop -c $inpcrd -r $min_rst -x $min_crd -o $min_out
#echo "Minimizacao finalizada com sucesso!"

#echo "Iniciando relaxamento parte 1:"
#echo " 300ps de pré-relaxamento NPT com restrição para proteína e ligante"
#$do_parallel -O -i $input_relax_part_1 -p $prmtop -c $min_rst -ref $min_rst -r $relax_rst_1 -x $relax_crd_1 -o $relax_out_1
#echo "Relaxamento parte 1 finalizada com sucesso!"

#echo "Iniciando relaxamento parte 2:"
#echo "300ps de relaxamento NPT com restrição para proteína"
#$do_cuda -O -i $input_relax_part_2 -p $prmtop -c $relax_rst_1 -ref $relax_rst_1 -r $relax_rst_2 -x $relax_crd_2 -o $relax_out_2
#echo "Relaxamento parte 2 finalizada com sucesso!"

#echo "Iniciando relaxamento parte 3:"
#echo "200ps de relaxamento NPT sem restrição para as cadeias laterais dos resíduos em até 5 Angstrons ao redor do ligante"
#$do_cuda -O -i $input_relax_part_3 -p $prmtop -c $relax_rst_2 -ref $relax_rst_2 -r $relax_rst_3 -x $relax_crd_3 -o $relax_out_3
#echo "Relaxamento parte 3 finalizada com sucesso!"

#echo "Iniciando relaxamento parte 4:"
#echo "200ps de relaxamento NPT sem restrição para os resíduos em até 5 Angstrons ao redor do ligante"
#$do_cuda -O -i $input_relax_part_4 -p $prmtop -c $relax_rst_3 -ref $relax_rst_3 -r $relax_rst_4 -x $relax_crd_4 -o $relax_out_4
#echo "Relaxamento parte 4 finalizada com sucesso!"

#echo "Iniciando pre-producao:"
#echo "5 ns"
#$do_cuda -O -i $input_pre_prod -p  $prmtop  -c $relax_rst_4 -ref $relax_rst_4  -r $pre_prod_rst -x $pre_prod_crd -o $pre_prod_out
#echo "Pre-producao finalizada com sucesso!"

echo "Iniciando producao:"
echo "100 ns"
$do_cuda -O -i $input_production -p  $prmtop  -c $pre_prod_rst -ref $pre_prod_rst  -r $production_rst -x $production_crd -o $production_out
echo "Producao finalizada com sucesso!"

#########################################################################################################################################

trap : 0

echo >&2 '
*************
***  FIM  *** 
*************
'
