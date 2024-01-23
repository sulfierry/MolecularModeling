#!/bin/bash

total=7 # Número total de iterações
echo "Iniciando o processo de remoção de água e íons..."

for i in {4..10}
do
    echo -n "[$i/$total] Processando pasta $i..."

    # Definindo os parâmetros
    parameters=$i/5cc8.prmtop
    trajetoria=$i/production.crd
    output=./water_remov/$i/5cc8_${i}_wr

    # Criando e escrevendo no arquivo de configuração do cpptraj
    echo "parm $parameters" > cpptraj.in
    echo "parm $parameters [top1]" >> cpptraj.in
    echo "trajin $trajetoria" >> cpptraj.in
    echo "autoimage" >> cpptraj.in
    echo "parmstrip :WAT,Cl-,Na+ parmindex 1" >> cpptraj.in
    echo "parmwrite out $output.prmtop parmindex 1" >> cpptraj.in
    echo "strip :WAT,Cl-,Na+" >> cpptraj.in
    echo "outtraj $output.pdb onlyframes 1" >> cpptraj.in
    echo "trajout $output.dcd" >> cpptraj.in
    echo "run" >> cpptraj.in

    # Executando cpptraj com o arquivo de configuração
    cpptraj < cpptraj.in

    # Removendo o arquivo de configuração
    rm cpptraj.in

    echo " Concluído."
done

echo "Processo concluído para todas as pastas."
