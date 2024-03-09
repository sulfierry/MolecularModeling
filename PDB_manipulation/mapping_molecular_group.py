

def map_to_molecular_group(atom_name, residue_name):
    residue_name = residue_name.upper()
    atom_name = atom_name.upper()

    # Verificação de Carbonos
    # Verificação de Carbonos
    if atom_name.startswith("C"):
        carbon_mappings = {
            ("VAL", ("CG1", "CG2")): "CR",
            ("LEU", ("CG", "CD1", "CD2")): "CR",
            ("ILE", ("CG1", "CG2", "CD1")): "CR",
            ("MET", ("CG", "SD", "CE")): "CR",
            ("PHE", ("CZ",)): "CB",
            ("TYR", ("CZ",)): "CB",
            ("TRP", ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2")): "CB",
            ("ASP", ("CG", "CD")): "CO2M",
            ("GLU", ("CG", "CD")): "CO2M",
            ("ARG", ("CG",)): "CGD"
        }
        for key, value in carbon_mappings.items():
            if residue_name == key[0] and atom_name in key[1]:
                return value
        if atom_name == "CA" and residue_name != "PRO":
            return "HC"
        
        else: return "CR"
        
    # Verificação de Oxigênios
    if atom_name.startswith("O"):
        if atom_name in ["OD1", "OE1", "OE2", "O"]:
            return "O=C"
        elif atom_name in ["OD2", "OG", "OG1", "OH"]:
            return "OR"
        elif residue_name in ["ASP", "GLU"] and atom_name in ["OD1", "OD2", "OE1", "OE2"]:
            return "O2CM"
        else: 
            return "OR"


    # Verificação de Nitrogênios
    if atom_name.startswith("N"):
        nitrogen_mappings = {
            ("LYS", "NZ"): "NR",
            ("ARG", "NE", "NH1", "NH2"): "NGD+",
            ("HIS", "ND1", "NE2"): "NPYD",
            ("ASN", "ND2"): "NC=O",
            ("GLN", "NE2"): "NC=O",
            ("TRP", "NE1"): "NPYL",
            ("PRO", "N"): "NC=C",
            ("CYS", "NSP"): "NSP"
        }
        for key, value in nitrogen_mappings.items():
            if residue_name == key[0]:
                if atom_name in key[1]:
                    return value
            else:
                return "NR"

    # Verificação de Hidrogênios
    if atom_name.startswith("H"):
        hydrogen_mappings = {
            ("SER", "HO"): "HOR",
            ("THR", "HO"): "HOR",
            ("LYS", "HNZ"): "HNR",
            ("ASP", "HOCO"): "HOCO",
            ("GLU", "HOCO"): "HOCO"
        }
        for key, value in hydrogen_mappings.items():
            if residue_name == key[0]:
                if atom_name in key[1]:
                    return value
        if atom_name == "HN" or atom_name[1:].isdigit():
            return "HNCO"
        else:
            return "H"

    # Verificação de Enxofres
    if atom_name.startswith("S"):
        sulfur_mappings = {
            ("MET", ["SD"]): "S",
            ("CYS", ["SG"]): "HS"
        }
        for key, value in sulfur_mappings.items():
            if residue_name == key[0]:
                if atom_name in key[1]:
                    return value
                
    if atom_name.startswith("MG"):
        return "MG+2"
                
    print(f"O átomo {atom_name} não foi reconhecido.")
    return "UNKNOWN"

