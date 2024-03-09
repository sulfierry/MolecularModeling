import math
import csv
from interaction import *


# Função para calcular a distância entre dois átomos
def calculate_distance(atom1, atom2):
    coord1 = atom1['coord']
    coord2 = atom2['coord']
    return math.sqrt(sum([(c1 - c2) ** 2 for c1, c2 in zip(coord1, coord2)]))

# Função para calcular o produto interno de dois vetores
def dot_product(v1, v2):
    return sum([a*b for a, b in zip(v1, v2)])

# Função para calcular a norma (magnitude) de um vetor
def norm(v):
    return math.sqrt(dot_product(v, v))

# Função para calcular o produto vetorial de dois vetores
def cross_product(v1, v2):
    return [
        v1[1]*v2[2] - v1[2]*v2[1],
        v1[2]*v2[0] - v1[0]*v2[2],
        v1[0]*v2[1] - v1[1]*v2[0]
    ]

# Função para calcular o ângulo entre três pontos
def calculate_angle(coord_A, coord_B, coord_C):

    # 1 - Obtemos os vetores BA e BC.
    BA = [a-b for a, b in zip(coord_A, coord_B)]
    BC = [c-b for c, b in zip(coord_C, coord_B)]

    # 2 - Calculamos o produto escalar (dot product) desses dois vetores.
    # 3 - Calculamos a magnitude (norma) de cada vetor.
    # 4 - O cosseno do ângulo entre os vetores é dado por: cos(θ) = BA⋅BC / |BA|⋅|BC|
    cosine_angle = dot_product(BA, BC) / (norm(BA) * norm(BC))
    
    # 5 - Usamos a função arco-cosseno (acos) para obter o ângulo θ em radianos.
    angle = math.acos(max(-1.0, min(1.0, cosine_angle)))

    # 6 - Convertemos o ângulo de radianos para graus.
    return math.degrees(angle)

# Função para calcular o diedro entre quatro pontos
def calculate_dihedral(coord_A, coord_B, coord_C, coord_D):
    
    # 1 - Calculamos os vetores BA, CB e DC.
    # 2 - O primeiro plano é definido pelos vetores BA e CB, e o segundo plano é definido pelos vetores CB e DC.

    BA = [a-b for a, b in zip(coord_A, coord_B)]
    CB = [b-c for b, c in zip(coord_B, coord_C)]
    DC = [c-d for c, d in zip(coord_C, coord_D)]

    # 3 - Calculamos os vetores normais a esses planos usando o produto vetorial (cross product). 
    # O vetor normal ao primeiro plano é:  N1 = BA x CB    
    # E o vetor normal ao segundo plano é: N2 = CB x DC.
    normal1 = cross_product(BA, CB)
    normal2 = cross_product(CB, DC)
    n1_norm = norm(normal1)
    n2_norm = norm(normal2)

    # 4 - Calculamos o cosseno do ângulo entre os vetores normais usando o produto escalar e as magnitudes:
    #        cos(θ) = N1⋅N2 / |N1|⋅|N2|
    # Além disso, para garantir que o diedro esteja no intervalo correto de -180° a 180°, 
    # levamos em consideração a direção do vetor formado pelo produto vetorial dos vetores normais e o vetor CB. 
    # Se a direção for negativa, invertemos o sinal do cosseno.
    if n1_norm != 0:
        normal1 = [n/n1_norm for n in normal1]
    if n2_norm != 0:
        normal2 = [n/n2_norm for n in normal2]
    cosine_angle = dot_product(normal1, normal2)
    direction = dot_product(cross_product(normal1, normal2), CB)
    if direction < 0:
        cosine_angle = -cosine_angle

    # 5 - Usamos a função arco-cosseno (acos) para obter o ângulo θ em radianos.
    # 6 - Convertemos o ângulo de radianos para graus.
    angle = math.acos(max(-1.0, min(1.0, cosine_angle)))

    return math.degrees(angle)



def calculate_angles_for_nearest_atoms(parsed_data, near_residues_dict):
    all_atoms = parsed_data['chains'] + parsed_data['cofactors'] + parsed_data['ligands']

    def find_two_nearest_atoms(target_atom, atoms_list):
        distances = [(atom, calculate_distance(target_atom, atom)) for atom in atoms_list if atom != target_atom]
        distances.sort(key=lambda x: x[1])
        return [distances[0][0], distances[1][0]]

    angles = []
    for entry in near_residues_dict:
        target_atom_serial = entry['molecule_atom_serial']
        target_atom = next(atom for atom in all_atoms if atom['serial_number'] == target_atom_serial)
        
        nearest_atoms = find_two_nearest_atoms(target_atom, all_atoms)
        
        angle = calculate_angle(nearest_atoms[0]['coord'], target_atom['coord'], nearest_atoms[1]['coord'])
        
        angles.append({
            'Target Atom'     : target_atom_serial,
            'Nearest Atoms'   : (nearest_atoms[0]['serial_number'], nearest_atoms[1]['serial_number']),
            'Angle (degrees)' : angle,
            'Atom Names'      : (nearest_atoms[0]['res_name'] + "("+ nearest_atoms[0]['name'] + ")" +  str(nearest_atoms[0]['res_seq']), 
                                target_atom['res_name']       + "("+ target_atom['name']      + ")" +  str(target_atom['res_seq']),
                                nearest_atoms[1]['res_name']  + "("+ nearest_atoms[1]['name'] + ")" +  str(nearest_atoms[1]['res_seq']))
        })


    return angles


def calculate_dihedrals_for_nearest_atoms(parsed_data, near_residues_dict):
    all_atoms = parsed_data['chains'] + parsed_data['cofactors'] + parsed_data['ligands']

    def find_three_nearest_atoms(target_atom, atoms_list):
        distances = [(atom, calculate_distance(target_atom, atom)) for atom in atoms_list if atom != target_atom]
        distances.sort(key=lambda x: x[1])
        return [distances[0][0], distances[1][0], distances[2][0]]

    dihedrals = []
    for entry in near_residues_dict:
        target_atom_serial = entry['molecule_atom_serial']
        target_atom = next(atom for atom in all_atoms if atom['serial_number'] == target_atom_serial)
        
        nearest_atoms = find_three_nearest_atoms(target_atom, all_atoms)
        
        dihedral = calculate_dihedral(nearest_atoms[0]['coord'], nearest_atoms[1]['coord'], nearest_atoms[2]['coord'], target_atom['coord'])
        
        dihedrals.append({
            'Target Atom': target_atom_serial,
            'Nearest Atoms': (nearest_atoms[0]['serial_number'], nearest_atoms[1]['serial_number'], nearest_atoms[2]['serial_number']),
            'Dihedral (degrees)': dihedral,
            'Atom Names': ("("+nearest_atoms[0]['name']+")" + str(nearest_atoms[0]['res_seq']), 
                           "("+nearest_atoms[1]['name']+")" + str(nearest_atoms[1]['res_seq']), 
                           "("+nearest_atoms[2]['name']+")" + str(nearest_atoms[2]['res_seq']),
                           "("+target_atom['name']+")"      + str(target_atom['res_seq']))
        })

    return dihedrals



def get_distance_between_atoms(parsed_data, atom1_idx, atom2_idx):
    atom1 = next(atom for atom in parsed_data['chains'] + parsed_data['cofactors'] + parsed_data['ligands'] if atom['serial_number'] == atom1_idx)
    atom2 = next(atom for atom in parsed_data['chains'] + parsed_data['cofactors'] + parsed_data['ligands'] if atom['serial_number'] == atom2_idx)
    return calculate_distance(atom1, atom2)



def set_output(near_residues_dict, ligand_residue_tuple, parsed_data, output_name):
    ligand_name, ligand_num, ligand_chain = ligand_residue_tuple
    
    # Ensure that angles_data and dihedrals_data have at least the same length as near_residues_dict
    angles_data = calculate_angles_for_nearest_atoms(parsed_data, near_residues_dict) + [{} for _ in range(len(near_residues_dict))]
    dihedrals_data = calculate_dihedrals_for_nearest_atoms(parsed_data, near_residues_dict) + [{} for _ in range(len(near_residues_dict))]
    
    interacting_molecules_count = 0

    with open(output_name, 'w', newline='') as file:
        writer = csv.writer(file)
        columns = ["Chain", "Nearby atoms", "Interaction", "Distance (Å)"]
        writer.writerow(columns)
        print("{:^5} {:^30} {:^20} {:^10}".format(*columns))

        for idx, entry in enumerate(near_residues_dict):
            chain_id = entry['chain']
            distance = entry['distance']
            molecule_atom = entry['molecule_atom']
            ligand_atom = entry['ligand_atom']
            aa_name = entry['molecule_name']
            aa_num = entry['molecule_number']

            atom1_str = f"{molecule_atom}({aa_name}{aa_num})"
            atom2_str = f"{ligand_atom}({ligand_name}{ligand_num})"
            nearby_atoms_str = f"{atom1_str:<12}-{atom2_str:>13}"

            probable_interaction = is_interaction(molecule_atom, ligand_atom, aa_name, distance)
            if probable_interaction not in ["Non-specific"]:
                interacting_molecules_count += 1

            angle_data = angles_data[idx]
            angle, angle_atoms, angle_names = round(angle_data.get('Angle (degrees)', ""), 2), angle_data.get('Nearest Atoms', ""), ", ".join(angle_data.get('Atom Names', ""))
            angle_distance = round(get_distance_between_atoms(parsed_data, *angle_atoms[:2]), 2) if 'Nearest Atoms' in angle_data else ""

            dihedral_data = dihedrals_data[idx]
            dihedral, dihedral_atoms, dihedral_names = round(dihedral_data.get('Dihedral (degrees)', ""), 2), dihedral_data.get('Nearest Atoms', ""), ", ".join(dihedral_data.get('Atom Names', ""))
            dihedral_distance = round(get_distance_between_atoms(parsed_data, *dihedral_atoms[1:3]), 2) if 'Nearest Atoms' in dihedral_data else ""

            # writer.writerow([chain_id, nearby_atoms_str, probable_interaction, round(distance, 2), round(angle,2), angle_distance, angle_names, dihedral, dihedral_distance, dihedral_names])
            # print("{:^5} {:^30} {:^20} {:^10.2f} {:^10} {:^15} {:^50} {:^15} {:^15} {:^50}".format(chain_id, nearby_atoms_str, probable_interaction, distance, angle, angle_distance, angle_names, dihedral, dihedral_distance, dihedral_names))

            writer.writerow([chain_id, nearby_atoms_str, probable_interaction, round(distance, 2)])
            print("{:^5} {:^30} {:^20} {:^10.2f}".format(chain_id, nearby_atoms_str, probable_interaction, distance))

    print("\nTotal number of interacting molecules:", interacting_molecules_count)
    print("\n")
    return f"Successfully processed and saved! Total interactions: {interacting_molecules_count}"
