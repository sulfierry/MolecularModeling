# Executar a conversão
obabel -i pdb 6kpf_top1.pdb -o mol2 -O 6kpf_top1.mol2

# Atribuir Tipos de Átomos:
antechamber -i 6kpf_top1.mol2 -fi mol2 -o output.mol2 -fo mol2 -c bcc -s 2

# Gerar Arquivo de Preparação (.frcmod - contém parâmetros adicionais necessários que podem não estar presentes no campo de força GAFF padrão)
parmchk2 -i output.mol2 -f mol2 -o ligand.frcmod

# tleap -f tleap.in
# source leaprc.gaff
# loadamberparams ligand.frcmod
# loadoff lig.lib 
# lig = loadmol2 output.mol2
# saveamberparm lig lig.prmtop lig.inpcrd
# savepdb lig output.pdb
# quit
