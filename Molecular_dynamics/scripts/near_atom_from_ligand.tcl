# Verifica e cria a pasta "vmd_frames" se ela não existir
if {![file exists vmd_frames]} {
    file mkdir vmd_frames
}

# Defina o nome do ligante
set ligand_name "LIG" ;# Substitua LIG pelo nome do resíduo do seu ligante

# Defina a distância de corte
set cutoff_distance 3.0

# Obtém o número total de quadros na trajetória
set num_frames [molinfo 0 get numframes]

# Itera sobre cada quadro da simulação
for {set frame_number 0} {$frame_number < 2000} {incr frame_number} {
    # Atualiza para o quadro atual
    molinfo 0 set frame $frame_number

    # Abre o arquivo para escrita para o quadro atual dentro da pasta "vmd_frames"
    set outfile [open "vmd_frames/detailed_interactions_frame_$frame_number.dat" w]

    # Cria uma seleção de todos os átomos do ligante
    set ligand_atoms [atomselect 0 "resname $ligand_name"]

    # Obtém informações detalhadas dos átomos do ligante
    foreach ligand_atom [$ligand_atoms get {index name}] {
        set ligand_index [lindex $ligand_atom 0]
        set ligand_atomname [lindex $ligand_atom 1]
        
        # Cria uma seleção de todos os átomos de cadeias laterais de aminoácidos próximos a este átomo do ligante
        set sidechain_atoms [atomselect 0 "sidechain within $cutoff_distance of index $ligand_index"]
        set sidechain_atoms_info [$sidechain_atoms get {index element resname resid name}]
        
        foreach sidechain_atom $sidechain_atoms_info {
            set sidechain_index [lindex $sidechain_atom 0]
            set sidechain_element [lindex $sidechain_atom 1]
            set sidechain_resname [lindex $sidechain_atom 2]
            set sidechain_resid [lindex $sidechain_atom 3]
            set sidechain_atomname [lindex $sidechain_atom 4]
            
            # Escreve as informações no arquivo, incluindo o nome do átomo do ligante
            puts $outfile "Frame $frame_number: Ligand atom $ligand_atomname index $ligand_index interacts with $sidechain_resname $sidechain_resid $sidechain_atomname"
        }
    }

    # Fecha o arquivo
    close $outfile
}

# Encerra o script VMD
quit
