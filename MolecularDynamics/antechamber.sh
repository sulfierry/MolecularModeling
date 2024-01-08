# Executar a conversão
obabel -i pdb 6kpf_top1.pdb -o mol2 -O 6kpf_top1.mol2

# Atribuir Tipos de Átomos:
antechamber -i 6kpf_top1.mol2 -fi mol2 -o output.mol2 -fo mol2 -c bcc -s 2

parmchk2 -i output.mol2 -f mol2 -o ligand.frcmod
