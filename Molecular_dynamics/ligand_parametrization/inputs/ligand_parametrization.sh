#!/bin/bash

# Converter PDB para MOL2 com AM1-BCC
antechamber -i ligand.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc

# O .frcmod contém parâmetros de força de campo adicionais que podem não estar presentes no GAFF (General Amber Force Field) padrão
parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod

python adjust_partial_charge.py

# Executar tleap com um arquivo de entrada tleap.in
tleap -f tleap.in
