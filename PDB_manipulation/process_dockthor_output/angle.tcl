set sele1 [atomselect top "resname LIG"]
set sele2 [atomselect top "resid 2 and name CA"]
set sele3 [atomselect top "resid 3 and name CA"]
set nf [molinfo top get numframes]
set outfile [open "angle.dat" a]
for {set i 0} {$i < $nf} {incr i} {

    $sele1 frame $i
    $sele2 frame $i
    $sele3 frame $i

    set com1 [measure center $sele1 weight mass]
    set com2 [measure center $sele2 weight mass]
    set com3 [measure center $sele3 weight mass]

    set avec [vecsub $com1 $com2]
    set amag [veclength $avec]

    set bvec [vecsub $com3 $com2]
    set bmag [veclength $bvec]

    set dotprod [vecdot $avec $bvec]

    set angle [expr 57.2958 * acos($dotprod / ($amag * $bmag))]

    puts $outfile "$i $angle"
}
close $outfile
quit
