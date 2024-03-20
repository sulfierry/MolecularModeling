set outfile [open "coordinates.dat" w]
set n [molinfo top get numframes]
for {set i 0} {$i < $n} {incr i} {
    set atomSel [atomselect top "resname LIG" frame $i] # ou uma seleção específica de átomos
    set coords [$atomSel get {x y z}]
    puts $outfile "$coords"
}
close $outfile
quit
