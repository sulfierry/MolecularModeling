#!/bin/bash

# Verifica se o número correto de argumentos foi fornecido
if [ "$#" -ne 3 ]; then
    echo "Uso: $0 input.prmtop input.crd nome_output"
    exit 1
fi

# Atribui os argumentos a variáveis
input_prmtop=$1
input_crd=$2
output_base_name=$3

total=7 # Número total de iterações
echo "Iniciando o processo de remoção de água e íons..."

# Criar a pasta "water_remov" se ela não existir
mkdir -p water_remov

for i in {4..10}
do
    echo -n "[$i/$total] Processando pasta $i..."

    # Define os caminhos de entrada e saída para cada iteração
    parameters="./$i/$input_prmtop"
    trajetoria="./$i/$input_crd"
    output="./water_remov/$i/${output_base_name}_${i}"

    # Criar subpasta dentro de "water_remov" se não existir
    mkdir -p ./water_remov/$i

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
