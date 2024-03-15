set outfile [open "distance.dat" a]
set n [molinfo top get numframes] 
set sum 0
for {set i 0} {$i < $n} {incr i} {
    set atom1 [atomselect top "resname LIG" frame $i] 
    set atom2 [atomselect top "resid 378" frame $i] 
    set com1 [measure center $atom1 weight mass]
    set com2 [measure center $atom2 weight mass]
    set dis [veclength [vecsub $com1 $com2]]
    set sum [expr $sum + $dis]
    puts $outfile "$i  $dis"
}
set Avg_dis [expr $sum/$n]
#puts $outfile "Avg Distance : $Avg_dis "
#puts "Avg Distance : $Avg_dis "
close $outfile
quit
