#!/bin/bash

total=10
echo "Iniciando o processo de remoção de água e íons..."

for i in {1..10}
do
    echo -n "[$i/$total] Processando pasta $i..."
    python script-remove-water-and-ions.py $i/5cc8.prmtop $i/production.crd ./water_remov/$i/5cc8_${i}_wr
    echo " Concluído."
done

echo "Processo concluído para todas as pastas."
