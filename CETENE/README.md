# Process Contextualization

## Calculation of Mulliken Charges
This step provides an initial charge distribution based on the electronic density obtained from quantum calculations. Mulliken charges may not be the most suitable for molecular dynamics simulations due to their dependence on the orbital basis used and other limitations.

## Application of the AM1-BCC Method
The AM1-BCC (Austin Model 1 - Bond Charge Correction) method is an improvement over Mulliken charges. It adjusts charges using a semi-empirical model to make them more appropriate for use in force fields like GAFF (General Amber Force Field). This method is known for providing more realistic charges for a wide range of organic molecules.

## Empirical Adjustment at the 3rd Decimal Place
The adjustment made in this step, though based on an empirical procedure, occurs after more rigorous calculations of the charges. The aim of this adjustment is to ensure that the total charge of the ligand is an integer (usually neutral or a specific ionic charge value). This adjustment is carried out in a way that minimizes the impact on the charges calculated by AM1-BCC.

# Ligand Parametrization Using Semi-Empirical Methods

This section outlines the process employed for the parametrization of ligands using semi-empirical methods, primarily within the AMBER software framework. The procedure involves several steps, each crucial for ensuring accurate representation of the ligand in molecular dynamics simulations.

## Steps Involved in the Parametrization Process

1. **Conversion of PDB to MOL2 Format**: 
   - Command: `antechamber -i ligand.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc`
   - Purpose: Convert the ligand file from PDB format to MOL2 format, suitable for further parametrization.

2. **Generation of Mulliken Charges using Antechamber**: 
   - Command: `antechamber -i ligand.mol2 -fi mol2 -fo ac -o ligand.ac -c mul`
   - Purpose: Calculate Mulliken charges for the ligand, providing an initial charge distribution.

3. **Application of AM1-BCC Method**: 
   - Command: `am1bcc -i ligand.ac -f ac -o ligand_bcc.ac -j 5`
   - Purpose: Employ the AM1-BCC method to refine the atomic charges based on a semi-empirical approach.

4. **Adjustment of Charges to Integer Values**: 
   - Command: `antechamber -i ligand_bcc.ac -fi ac -o ligand_charges_values.mol2 -fo mol2 -nc <desired_charge>`
   - Purpose: Adjust the total charge of the ligand to an integer value, ensuring charge neutrality or a specific ionic state.

5. **Empirical Charge Distribution Among Hydrogen Atoms**: 
   - Script: `verify_partial_charge.py`
   - Purpose: Empirically distribute charges among hydrogen atoms to ensure the final total charge of the ligand is an integer. This script fine-tunes the partial charges of hydrogen atoms to align with the molecule's total charge.

6. **Generation of Additional Force Field Parameters (.frcmod File)**: 
   - Command: `parmchk2 -i ligand_charges_values.mol2 -f mol2 -o ligand.frcmod`
   - Purpose: Generate a .frcmod file containing additional force field parameters that might not be present in the standard GAFF (General Amber Force Field).

7. **Preparation of Ligand Library and Complex for Simulation**: 
   - Commands: 
     - `source leaprc.gaff`
     - `LIG = loadmol2 ligand_charges_values.mol2`
     - `saveoff LIG ligand.lib`
     - `saveamberparm LIG ligand.prmtop ligand.rst7`
   - Purpose: Prepare the ligand library and parameter files necessary for the simulation, using the tleap tool of AMBER.

8. **Assembly of the Protein-Ligand Complex**: 
   - Commands: 
     - `source leaprc.protein.ff19SB`
     - `source leaprc.water.opc`
     - `loadamberparams ligand.frcmod`
     - `loadoff ligand.lib`
     - `complex = loadpdb complex.pdb`
     - `solvateOct complex OPCBOX 15.0`
     - `addIons complex Cl- 14`
     - `saveamberparm complex complex.prmtop complex.rst7`
   - Purpose: Construct the full protein-ligand simulation system, including solvation and ion placement, and
