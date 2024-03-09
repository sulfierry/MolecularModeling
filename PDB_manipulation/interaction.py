from mmff94_molecular_group import *
from mmff94_vdw import *
from molecule_class import *
from mapping_molecular_group import *


# Distance threshold for hydrophobic interactions
hydrophobic_distance_threshold = 4.0

# Define a function to check for potential hydrogen bonds
def is_interaction(atom1_name, atom2_name, residue_name, distance):

    # Check for hydrophobic interactions
    if residue_name in hydrophobic_residues:
        if atom1_name in hydrophobic_atoms or atom2_name in hydrophobic_atoms:
            if distance <= hydrophobic_distance_threshold:
                return "Hydrophobic"
            
        # Specifically check for a carbon-hydrogen interaction
        if (atom1_name.startswith("C") and atom2_name == "H") or (atom2_name.startswith("C") and atom1_name == "H"):
            if distance <= hydrophobic_distance_threshold:
                return "Hydrophobic"
            
    # Check for ionic interaction
    if distance < 4:
        if residue_name in ionic_interactions:
            if atom1_name.startswith(tuple(ionic_interactions.get(residue_name, []))) or \
                atom2_name.startswith(tuple(ionic_interactions.get(residue_name, []))):
                return "Ionic interaction"
            
        elif (atom1_name in ionic_atoms['positive'] and atom2_name in ionic_atoms['negative']) or \
           (atom2_name in ionic_atoms['positive'] and atom1_name in ionic_atoms['negative']):
            return "Ionic interaction"
        
    # Check for hydrogen bond      
    if distance < 3.7:   
        if atom1_name.startswith(tuple(hydrogen_bond_acceptors)) and \
           atom2_name.startswith(tuple(hydrogen_bond_acceptors)) and \
           atom1_name != atom2_name:        
            return "Hydrogen bond"

    # Assuming lennard_jones_potential is defined elsewhere in your code
    v_lj = lennard_jones_potential(atom1_name, atom2_name, residue_name, distance)

    # Check for potential van der Waals interaction based on Lennard-Jones potential
    # Valor Conservador: Se você deseja ser mais conservador e focar apenas nas interações 
    # mais fortes de van der Waals, pode considerar um valor limiar de −0.5 kcal/mol ou mais negativo.
    
    # Valor Moderado: Um valor de −0.2 a −0.3 kcal/mol pode ser uma abordagem intermediária, 
    # onde você identifica interações que têm uma contribuição notável, mas não são extremamente fracas.

    # Análise Detalhada: Se o objetivo é uma análise mais detalhada e abrangente das interações, 
    # incluindo as mais fracas, então −0.1 kcal/mol ou até um pouco mais positivo pode ser aceitável. 
    # No entanto, essas interações devem ser interpretadas com cautela e corroboradas com outras evidências ou análises.

    if v_lj < -0.1:
        return "van der Waals"

    return "Non-specific"


def lennard_jones_potential(atom1_name, atom2_name, residue, r):

    """Calculate Lennard-Jones potential between two atoms based on MMFF types."""
    # Get MMFF types for the atoms
    mapped_group1 = map_to_molecular_group(atom1_name, residue)
    mapped_group2 = map_to_molecular_group(atom2_name, residue)

    if not mapped_group1 or not mapped_group2: 
        return 0  # or handle this case as required

    type1 = mmff94_molecular_group[mapped_group1]['PRIMARY MMF TYPE']
    type2 = mmff94_molecular_group[mapped_group2]['PRIMARY MMF TYPE']

    # Get epsilon and sigma values for the atoms
    epsilon1 = float(mmff94_vdw[type1]['alpha-i'])
    epsilon2 = float(mmff94_vdw[type2]['alpha-i'])

    sigma1 = float(mmff94_vdw[type1]['N-i'])
    sigma2 = float(mmff94_vdw[type2]['N-i'])

    # Combine the epsilon and sigma values
    # Este cálculo refere-se à combinação dos parâmetros de profundidade do poço de energia epsilon para dois átomos 
    # diferentes quando se modela uma interação via potencial de Lennard-Jones. 
    # A combinação geométrica (média geométrica) é comum para este parâmetro.
    epsilon_combined = (epsilon1 * epsilon2) ** 0.5

    # Este cálculo refere-se à combinação dos parâmetros de distância de sigma para os mesmos dois átomos. 
    # Sigma é geralmente interpretado como a distância em que o potencial interatômico entre dois átomos neutros é zero.
    # A combinação aritmética (média aritmética) é típica para este parâmetro.
    sigma_combined = (sigma1 + sigma2) / 2.0

    # Calculate the Lennard-Jones potential
    return 4 * epsilon_combined * ((sigma_combined / r)**12 - (sigma_combined / r)**6)
