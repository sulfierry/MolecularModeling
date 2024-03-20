set outfile [open "coordinates.dat" a]
set n [molinfo top get numframes]
for {set i 0} {$i < $n} {incr i} {
    set atomSel [atomselect top "resname LIG" frame $i]
    # Agora, também obtemos os nomes dos átomos
    set coordsList [$atomSel get {x y z}]
    set atomNames [$atomSel get name]
    foreach coords $coordsList name $atomNames {
        # coords é uma lista de três elementos: x, y, z
        # name é o nome do átomo
        # Formata a saída como "x y z nome" e escreve no arquivo
        puts $outfile "[join $coords " "] $name"
    }
}
close $outfile
quit
