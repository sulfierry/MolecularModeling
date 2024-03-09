from geometric_manipulation import *
from mmff94_molecular_group import *
from molecule_class import *
from interaction import *
from mmff94_vdw import *
from parse_pdb import *



if __name__ == "__main__":

    # Arquivos de entrada e saida a serem fornecidos
    input_pdb      = "./3c9t.pdb"  #sys.argv[1]    # EXAMPLE.pdb
    input_molecule = "ACP"         #sys.argv[2]    # ATP
    output_name    = "out.csv"     #sys.argv[3]    # ATP_OUT (csv)
 
    # Distance from the selected molecule
    treshold_distance = 4.0

    # executa e salva o resultados para a classificacao dos contatos
    input_pdb      = parse_pdb(input_pdb)
    ligand_residue = find_molecule(input_pdb, input_molecule)
    near_residues  = verify_near_residues(input_pdb, ligand_residue, treshold_distance)
    set_output(near_residues, ligand_residue, input_pdb, output_name)
