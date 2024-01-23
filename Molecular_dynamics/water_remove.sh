#!/bin/bash

# Uso: ./script.sh xxx.prmtop nome-da-trajetoria.mdcrd nome-para-os-outputs

parameters=$1
trajetoria=$2
output=$3

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
