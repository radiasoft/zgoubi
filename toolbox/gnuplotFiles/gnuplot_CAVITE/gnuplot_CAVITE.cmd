
reset
 
set xlabel "time [ms]" font "roman,18"
set x2label "turn numb. [{/Symbol \264}10^3]" font "roman,18"
set ylabel "RF voltage [MV]" font "roman,18"
set y2label "E_s [MeV]" font "roman,18"

set title "RF voltage and particle energy" font "roman,20"
 
set xtics nomirror font "roman,11"
set x2tics nomirror font "roman,11"
set ytics nomirror font "roman,11"
set y2tics nomirror font "roman,11"
 
set key  t c font "roman,10"  

c = 2.99792458e8
Circ = 755.8699
Trev = Circ / c
V2MV = 1e-3
s2ms = 1e3

t1 = 0 ; t2 = 9. 
set xrange [t1:t2]
set x2range [t1/Trev/s2ms /1000.:t2/Trev/s2ms/1000.]

plot \
    'zgoubi.CAVITE.Out'  u ($9 *Trev *s2ms):($21) axes x1y1 w l tit "V" ,\
    'zgoubi.CAVITE.Out'  u ($9 *Trev *s2ms):($17 *V2MV) axes x1y2 w l tit "E_s" 

set terminal postscript eps blacktext color enh size 8cm,5cm "Times-Roman" 12 
 set output "gnuplot_CAVITE.eps" 
 replot 
 set terminal X11 
 unset output 

pause 2
exit


