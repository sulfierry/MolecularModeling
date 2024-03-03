# vmd -dispdev text -e script.tcl 
# Find the number of atoms within 2-3 A of the ligand. 

# Define o nome do ligante
set ligand_name "LIG"

# Cria uma seleção de todos os resíduos de proteína próximos ao ligante
set closer_residue [atomselect 0 "protein within 3.00 of resname $ligand_name"]

# Abre o arquivo para escrita
set outfile [open "residues_near_ligand.dat" w]

# Obtém o número total de quadros na trajetória
set num_frames [molinfo 0 get numframes]

# Itera sobre cada quadro da simulação
for {set frame_number 0} {$frame_number < $num_frames} {incr frame_number} {
    # Atualiza a seleção para o quadro atual
    $closer_residue frame $frame_number
    $closer_residue update

    # Obtém os resíduos únicos próximos ao ligante neste quadro
    set residues [$closer_residue get {resname resid}]
    set unique_residues [lsort -unique $residues]

    # Escreve o número do quadro e o número de contatos únicos
    puts $outfile "Frame $frame_number: [llength $unique_residues] unique residues near $ligand_name"

    # Para cada resíduo único, imprime sua identificação
    foreach {resname resid} $unique_residues {
        puts $outfile "    Residue $resname $resid"
    }
}

# Fecha o arquivo
close $outfile

# Encerra o script VMD
quit
