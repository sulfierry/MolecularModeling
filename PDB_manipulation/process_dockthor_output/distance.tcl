set outfile [open "distance.dat" a]
set n [molinfo top get numframes]
set sum 0
for {set i 0} {$i < $n} {incr i} {
    set atom1 [atomselect top "index 9736" frame $i]
    set atom2 [atomselect top "resid 486" frame $i]
    set com1 [measure center $atom1]
    set com2 [measure center $atom2]
    set dis [veclength [vecsub $com1 $com2]]
    set sum [expr $sum + $dis]
    puts $outfile "$i  $dis"
}
close $outfile
quit
