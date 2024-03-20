set outfile [open "coordinates.dat" a]
set n [molinfo top get numframes]
for {set i 0} {$i < $n} {incr i} {
    set atomSel [atomselect top "resname LIG" frame $i]
    # Obtem as coordenadas, os nomes, os índices dos átomos e o índice do frame atual
    set coordsList [$atomSel get {x y z}]
    set atomNames [$atomSel get name]
    set atomIndices [$atomSel get index]
    foreach coords $coordsList name $atomNames index $atomIndices {
        # coords é uma lista de três elementos: x, y, z
        # name é o nome do átomo, e index é o índice do átomo
        # i é o índice do frame atual
        # Formata a saída como "frame x y z nome índice" e escreve no arquivo
        puts $outfile "$i [join $coords " "] $name $index"
    }
}
close $outfile
quit
