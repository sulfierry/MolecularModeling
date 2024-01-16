# 1. Converter PDB para MOL2
antechamber -i 6pt0_lig.pdb -fi pdb -o input.mol2 -fo mol2 -c bcc

# 2. Gerar Cargas Mulliken com o antechamber
antechamber -i input.mol2 -fi mol2 -fo ac -o output.ac -c mul

# 3. Aplicar o AM1-BCC com o am1bcc
am1bcc -i output.ac -f ac -o updated_output.ac -j 5

# 4. Ajustar Cargas para um Número Inteiro <desired_charge = 0 [neste exemplo]> 
antechamber -i updated_output.ac -fi ac -o ligand.mol2 -fo mol2 -nc <desired_charge>

# 5. Gerar .frcmod (contém parâmetros adicionais necessários que podem não estar presentes no campo de força GAFF padrão)
parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod

# 6. Gerar biblioteca do ligante (tleap)
source leaprc.gaff
LIG = loadmol2 ligand_adjusted.mol2
saveoff LIG ligand.lib
saveamberparm SUS sustiva.prmtop sustiva.rst7

# 7. Gerar complexo (tleap)
source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.gaff
loadamberparams ligand.frcmod
loadoff ligand.lib
complex = loadpdb 6pt0_complex.pdb
charge complex
addIons complex Cl- 14
solvateOct complex OPCBOX 15.0
check complex
SaveAmberParm complex complex.prmtop complex.rst7
