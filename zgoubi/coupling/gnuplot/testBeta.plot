set terminal postscript eps blacktext color enh size 15cm,7cm "Times-Roman" 17
set output "test-betaI.eps"

set title"GENERALIZED TWISS' PARAMETER {/Symbol b} FOR THE EIGENMODE I"
set xlabel"Arclength (in m)"
set xrange[0:200]
set ylabel"{/Symbol b}_I (in m)"

set key spacing 3
set ytics nomirror
set grid
set key outside

plot [][] \
 "madx.out" u ($1):($2) w lp pt 1 ps 1. tit "{/Symbol b}@^{/Times-Roman= 10 MADX}_I", \
 "twiss.out" u ($3):($4) w lp pt 1 ps 1. tit "{/Symbol b}@^{/Times-Roman= 10 ET}_I"#, \


######################################################################


set terminal postscript eps blacktext color enh size 15cm,7cm "Times-Roman" 17
set output "test-betaII.eps"

set title"GENERALIZED TWISS' PARAMETER {/Symbol b} FOR THE EIGENMODE II"
set xlabel"Arclength (in m)"
set xrange[0:200]
set ylabel"{/Symbol b}_{II} (in m)"

set key spacing 3
set ytics nomirror
set grid
set key outside

plot [][] \
 "madx.out" u ($1):($3) w lp pt 2 ps 1. tit "{/Symbol b}@^{/Times-Roman= 10 MADX}_{II}", \
 "twiss.out" u ($3):($5) w lp pt 2 ps 1. tit "{/Symbol b}@^{/Times-Roman= 10 ET}_{II}"#, \


exit
