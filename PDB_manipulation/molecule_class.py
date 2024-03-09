
# Amino acids and nucleic acid bases classification
molecule_class = {

    # Amino Acids
    "ALA": "hydrophobic",
    "ILE": "hydrophobic",
    "LEU": "hydrophobic",
    "VAL": "hydrophobic",
    "PHE": "hydrophobic",
    "PRO": "hydrophobic",
    "TRP": "hydrophobic",
    "MET": "hydrophobic",
    "GLY": "hydrophobic",
    "CYS": "polar charge: 0",
    "SER": "polar charge: 0",
    "THR": "polar charge: 0",
    "TYR": "polar charge: 0",
    "ASN": "polar charge: 0",
    "GLN": "polar charge: 0",
    "HIS": "polar charge: +",
    "LYS": "polar charge: +",
    "ARG": "polar charge: +",
    "ASP": "polar charge: -",
    "GLU": "polar charge: -",

    # Nucleic Acid Bases
    "A": "adenine base",
    "C": "cytosine base",
    "G": "guanine base",
    "T": "thymine base",   # DNA
    "U": "uracil base",    # RNA

    # Cofactors
    "MG" : "cofactor",
    "K"  : "cofactor",
    "CA" : "cofactor",
    "NA" : "cofactor",
    "HOH": "partial polar charge",

    # ligands
    "TPS": "ligand",
    "ACP": "ligand",
    "TMP": "ligand",
    "TPP": "ligand",
    "ANP": "ligand",   
}

ionic_interactions = {

        "NA": ["F", "Cl", "Br", "I", "O"],
        "MG": ["F", "Cl", "Br", "I", "O"],
        "K" : ["F", "Cl", "Br", "I", "O"],
        "CA": ["F", "Cl", "Br", "I", "O"],
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

