import sys
import math


# Parser for PDB structure with consideration for protein structure, ligands, and cofactors
def parse_pdb(pdb_file):
    """
    Parses a PDB file to extract chains, cofactors, ligands, and atom details.

    Parameters:
        pdb_file (str): Path to the PDB file.

    Returns:
        dict: A dictionary containing lists:
              - 'chains': List of ATOM details.
              - 'cofactors': List of HETATM details matching cofactor names.
              - 'ligands': List of HETATM details matching ligand names.
    """
    
    def extract_atom_details(line):
        """Helper function to extract atom details from a line."""
        return {
            'serial_number': int(line[6:11].strip()),
            'name': line[12:16].strip(),
            'alt_loc': line[16].strip(),
            'res_name': line[17:20].strip(),
            'chain_id': line[21].strip(),
            'res_seq': int(line[22:26].strip()),
            'icode': line[26].strip(),
            'coord': [
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ],
            'occupancy': float(line[54:60].strip()),
            'temp_factor': float(line[60:66].strip()),
            'element': line[76:78].strip(),
            'charge': line[78:80].strip()
        }

    chains = []
    cofactors = []
    ligands = []

    cofactor_names = ["MG", "ZN", "CA", "K", "NA", "FE", "CL", "HOH"]
    ligand_names = ["ACP", "TPS", "TMP", "TPP", "LIG"]

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chains.append(extract_atom_details(line))
            elif line.startswith('HETATM'):
                residue_name = line[17:20].strip()
                details = extract_atom_details(line)
                if residue_name in cofactor_names:
                    cofactors.append(details)
                elif residue_name in ligand_names:
                    ligands.append(details)

    return {
        'chains': chains,
        'cofactors': cofactors,
        'ligands': ligands
    }


def format_line(atom_data, atom_type="ATOM  "):
    return (
        f"{atom_type:6s}{atom_data['serial_number']:5d} {atom_data['name']:<4s} {atom_data['alt_loc']:1s}{atom_data['res_name']:<3s} "
        f"{atom_data['chain_id']:1s}{atom_data['res_seq']:4d}{atom_data['icode']:1s}   "
        f"{atom_data['coord'][0]:8.3f}{atom_data['coord'][1]:8.3f}{atom_data['coord'][2]:8.3f}"
        f"{atom_data['occupancy']:6.2f}{atom_data['temp_factor']:6.2f}          "
        f"{atom_data['element']:^2s}{atom_data['charge']:2s}\n"
    )

def print_pdb_structure(pdb_dict):
    """
    Prints the PDB structure from the dictionary in an organized manner.

    Parameters:
        pdb_dict (dict): Dictionary containing the parsed PDB lists.

    Returns:
        None: Simply prints the PDB structured data.
    """
    # Print the chains first
    for atom in pdb_dict['chains']:
        print(format_line(atom), end='')

    # Print cofactors
    print("TER")
    for atom in pdb_dict['cofactors']:
        print(format_line(atom, "HETATM"), end='')

    # Print ligands
    print("TER")
    for atom in pdb_dict['ligands']:
        print(format_line(atom, "HETATM"), end='')

    print("END")


def find_molecule(pdb_dict, molecule_name):
    """
    Searches for a molecule in the parsed PDB dictionary based on its name.

    Parameters:
        pdb_dict (dict): Dictionary containing parsed PDB data.
        molecule_name (str): Name of the molecule to search for.

    Returns:
        tuple: (Name of the molecule, residue sequence number) or None if not found.
    """

    chain_select = None

    # Check if the chain_select was provided as an argument
    if len(sys.argv) > 4:
        chain_select = str(sys.argv[4])

    # Combine all atom lists
    all_atoms = pdb_dict['chains'] + pdb_dict['cofactors'] + pdb_dict['ligands']

    # Check if any atom with the desired molecule name is present
    for atom in all_atoms:
        if atom['res_name'] == molecule_name:
            if chain_select is None or atom['chain_id'] == chain_select:
                return (atom['res_name'], atom['res_seq'], atom['chain_id'])
            else:
                atom['chain_id'] = chain_select
                return (atom['res_name'], atom['res_seq'], atom['chain_id'])

    return None


def verify_near_residues(input_pdb, ligand_residue, treshold_distance):
    ligand_atoms = [atom for atom in input_pdb['ligands'] 
                    if (atom['res_name'], atom['res_seq'], atom['chain_id']) == ligand_residue]

    # Combine all atom lists
    all_atoms = input_pdb['chains'] + input_pdb['cofactors'] + input_pdb['ligands']

    # Filter out the ligand atoms from all_atoms
    all_atoms = [atom for atom in all_atoms if (atom['res_name'], atom['res_seq'], atom['chain_id']) != ligand_residue]

    near_residues_dict = []

    for atom in all_atoms:
        # Find the closest ligand atom to the current atom
        min_distance, closest_ligand_atom = min((calculate_distance(atom, ligand_atom), ligand_atom) for ligand_atom in ligand_atoms)

        if min_distance <= treshold_distance:
            info = {
                'molecule_name': atom['res_name'],
                'molecule_number': atom['res_seq'],
                'chain': atom['chain_id'],
                'distance': min_distance,
                'molecule_atom': atom['name'],
                'ligand_atom': closest_ligand_atom['name'],
                'molecule_atom_serial': atom['serial_number'],  # Added the serial number for molecule atom
                'ligand_atom_serial': closest_ligand_atom['serial_number']  # Added the serial number for ligand atom
            }
            near_residues_dict.append(info)

    # Sort the list based on distance
    near_residues_dict.sort(key=lambda x: x['distance'])

    return near_residues_dict

def calculate_distance(atom1, atom2):
    coord1 = atom1['coord']
    coord2 = atom2['coord']
    return math.sqrt(sum([(c1 - c2) ** 2 for c1, c2 in zip(coord1, coord2)]))
