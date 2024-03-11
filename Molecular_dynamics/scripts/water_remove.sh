#!/bin/bash

# Verifica se o número correto de argumentos foi fornecido
if [ "$#" -ne 4 ]; then
    echo "Uso: $0 input_name.prmtop input_name.crd output_name num_pastas"
    echo "Example: ./water_remov.sh abc.prmtop abc.crd abc_wr 10"
    exit 1
fi

# Atribui os argumentos a variáveis
input_prmtop=$1
input_crd=$2
output_base_name=$3
num_pastas=$4

echo "Iniciando o processo de remoção de água e íons em $num_pastas pastas..."

# Criar a pasta "water_remov" se ela não existir
mkdir -p water_remov

process() {
    i=$1
    input_prmtop=$2
    input_crd=$3
    output_base_name=$4

    echo -n "[$i] Processando pasta $i..."

    # Define os caminhos de entrada e saída para cada iteração
    parameters="./$i/$input_prmtop"
    trajetoria="./$i/$input_crd"
    output="./water_remov/$i/${output_base_name}_${i}"

    # Criar subpasta dentro de "water_remov" se não existir
    mkdir -p "./water_remov/$i"

    # Criando e escrevendo no arquivo de configuração do cpptraj
    {
        echo "parm $parameters"
        echo "parm $parameters [top1]"
        echo "trajin $trajetoria"
        echo "autoimage"
        echo "parmstrip :WAT,Cl-,Na+ parmindex 1"
        echo "parmwrite out $output.prmtop parmindex 1"
        echo "strip :WAT,Cl-,Na+"
        echo "outtraj $output.pdb onlyframes 1"
        echo "trajout $output.dcd"
        echo "run"
    } > cpptraj.in

    # Executando cpptraj com o arquivo de configuração
    cpptraj < cpptraj.in

    # Removendo o arquivo de configuração
    rm cpptraj.in

    echo " Concluído."
}

export -f process

# Verifica se o comando parallel está disponível
if command -v parallel >/dev/null 2>&1; then
    # Usando parallel para executar a função process em paralelo
    parallel process ::: $(seq 1 $num_pastas) ::: $input_prmtop ::: $input_crd ::: $output_base_name
else
    # Executa de maneira serial se parallel não estiver disponível
    for i in $(seq 1 $num_pastas); do
        process $i $input_prmtop $input_crd $output_base_name
    done
fi

echo "Processo concluído para todas as $num_pastas pastas."
