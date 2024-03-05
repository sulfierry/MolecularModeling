from Bio.PDB import PDBParser, NeighborSearch
import sys
import csv


van_der_waals_radii = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 
    'F': 1.47, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98, 
    'Ca': 2.31, 'Fe': 2.01, 'Cu': 1.40, 'Zn': 1.39,
    'P': 1.80, 'S': 1.80, 'Se': 1.90, 
    # Adicionando mais halogênios
    'He': 1.40, 'Ne': 1.54, 'Ar': 1.88, 'Kr': 2.02, 'Xe': 2.16, 'Rn': 2.20, 
    # Outros átomos comuns
    'Na': 2.27, 'Mg': 1.73, 'K': 2.75,   # Valores exemplares para sódio, magnésio, potássio, ferro
    # Adicionando carbono alfa (CA) como um caso especial, geralmente tratado como carbono
    'CA': 1.70, # Usando o mesmo valor que para o carbono, mas pode precisar de ajuste dependendo do contexto
}


hydrogen_bond_acceptors = [
    "O",
    "N",
    "F",
    "H"
    ]

# Hydrophobic residues and atoms for identifying hydrophobic interactions
hydrophobic_residues = [
    "ALA",
    "VAL",
    "ILE",
    "LEU",
    "MET",
    "PHE",
    "TRP",
    "PRO",
    "TYR"
    ]


hydrophobic_atoms = [
    
    # Alanina
    "CB",
    # Valina
    "CB", "CG1", "CG2",
    # Isoleucina
    "CB", "CG1", "CG2", "CD1",
    # Leucina
    "CB", "CG", "CD1", "CD2",
    # Metionina
    "CB", "CG", "SD", "CE",
    # Prolina
    "CB", "CG", "CD",
    # Fenilalanina
    "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ",
    # Triptofano
    "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2",
    # Tirosina
    "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ",
    "H",
    # Átomos de carbono em grupos alquila e anéis alifáticos
    "C",
    # Nome comum para carbonos em anéis aromáticos
    "CA",
    # Outros nomes para carbonos em diferentes ambientes
    "CH",
    # Carbonos numerados, comuns em moléculas pequenas
    "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9",
    # Flúor
    "F",
    # Cloro
    "Cl",
    # Bromo
    "Br",
    # Iodo
    "I"
]

ionic_atoms = {
    'positive': ['NH1', 'NH2', 'NZ', 'ND1', 'NE2'],  # Átomos carregados positivamente
    'negative': ['OD1', 'OD2', 'OE1', 'OE2']  # Átomos carregados negativamente
}


# Distance threshold for hydrophobic interactions
hydrophobic_distance_threshold = 4.0

def is_hydrogen_bond(donor_atom, acceptor_atom, distance):
    # Verifica se um dos átomos é hidrogênio e se está ligado a um doador de hidrogênio (N ou O)
    if distance <= 3.5:
        if (donor_atom.element == 'H' and acceptor_atom.element in ['N', 'O']) or \
           (acceptor_atom.element == 'H' and donor_atom.element in ['N', 'O']):
            # Adicionalmente, verificar se o átomo de hidrogênio está ligado a um doador (N ou O)
            # Isso pode requerer acesso à estrutura molecular para verificar ligações
            return True
    return False

def is_hydrophobic_interaction(atom1, atom2, distance):
    # Critérios simplificados para interações hidrofóbicas
    hydrophobic_elements = ['C', 'S']
    if atom1.element in hydrophobic_elements and atom2.element in hydrophobic_elements:
        if distance <= 5.0:  # Distância em Angstroms
            return True
    return False

def is_van_der_waals_interaction(atom1, atom2, distance):
    # Obtém os raios de van der Waals para cada átomo
    radius1 = van_der_waals_radii.get(atom1.element, 0)
    radius2 = van_der_waals_radii.get(atom2.element, 0)

    # Calcula a distância máxima para uma interação de van der Waals
    max_distance = radius1 + radius2

    # Verifica se a distância real está dentro do limite
    return distance <= max_distance


def find_interactions(structure, ligand_atoms, threshold):
    ns = NeighborSearch(list(structure.get_atoms()))
    interactions = []
    for ligand_atom in ligand_atoms:
        for nearby_atom in ns.search(ligand_atom.coord, threshold, level='A'):
            if nearby_atom not in ligand_atoms:
                distance = ligand_atom - nearby_atom
                interaction_type = 'Unknown'

                # Ajuste na detecção de pontes de hidrogênio
                if ((ligand_atom.element == 'H' and nearby_atom.element in ['N', 'O']) or 
                    (nearby_atom.element == 'H' and ligand_atom.element in ['N', 'O'])) and distance <= 3.5:
                    interaction_type = 'Hydrogen bond'
                # Refinamento na detecção de interações hidrofóbicas
                elif ligand_atom.get_parent().get_resname() in hydrophobic_residues and nearby_atom.get_parent().get_resname() in hydrophobic_residues:
                    if distance <= hydrophobic_distance_threshold:
                        interaction_type = 'Hydrophobic'
                # Interações de van der Waals como fallback
                # Verifica interações de van der Waals
                if is_van_der_waals_interaction(ligand_atom, nearby_atom, distance):
                    interaction_type = 'van der Waals'

                interactions.append({
                    'ligand_atom': ligand_atom.get_name(),
                    'residue': nearby_atom.get_parent().get_resname(),
                    'residue_number': nearby_atom.get_parent().id[1],
                    'target_atom': nearby_atom.get_name(),
                    'distance': distance,
                    'interaction_type': interaction_type
                })
    return interactions

def save_to_csv(filename, interactions):
    with open(filename, 'w', newline='') as csvfile:
        fieldnames = ['residue_number', 'residue', 'ligand_atom', 'target_atom', 'distance', 'interaction_type']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for interaction in interactions:
            writer.writerow(interaction)

def print_structured_output(interactions):
    print("{:^12} {:^12} {:^15} {:^15} {:^20} {:^10}".format("Residue Num", "Residue", "Residue Atom", "Ligand Atom", "Interaction Type", "Distance (Å)"))
    print("-" * 85)
    for interaction in interactions:
        print("{:^12} {:^12} {:^15} {:^15} {:^20} {:^10.2f}".format(
            interaction['residue_number'],
            interaction['residue'],
            interaction['target_atom'],
            interaction['ligand_atom'],
            interaction['interaction_type'],
            interaction['distance']
        ))

def main(pdb_filename, ligand_resname, ligand_resnum, output_filename):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_filename)
    ligand_atoms = [atom for atom in structure.get_atoms() if atom.get_parent().get_resname() == ligand_resname and atom.get_parent().id[1] == ligand_resnum]
    interactions = find_interactions(structure, ligand_atoms, 5.0)  # Threshold distance set to 5.0 Angstroms
    
    # Aqui, você pode escolher entre salvar em um arquivo CSV ou imprimir no prompt.
    # Para salvar, você pode descomentar a próxima linha:
    # save_to_csv(output_filename, interactions)
    
    # Para imprimir no prompt de forma estruturada:
    print_structured_output(interactions)

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python script.py <PDB_FILE> <LIGAND_RESNAME> <LIGAND_RESNUM> <OUTPUT_CSV>")
    else:
        pdb_file = sys.argv[1]
        ligand_name = sys.argv[2]
        ligand_number = int(sys.argv[3])  # Ensure ligand number is an integer
        output_file = sys.argv[4]
        main(pdb_file, ligand_name, ligand_number, output_file)
