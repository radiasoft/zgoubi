
 set title "Orbits, from zgoubi.OPTICS.out"    font "roman,16"   # offset 0,+.7    

 set xlabel "s [m]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "x, y [m]"             font "roman,13"   #offset -0,-1 rotate by -20 
# set y2label "{/Symbol h}_x, {/Symbol h}_y"          font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12"   mirror
# set x2tics  font "roman,12" mirror
 set ytics  font "roman,12"   mirror      #offset 0,-.6
# set y2tics  font "roman,12" nomirror      #offset 0,-.6

set key b l maxrows 1 width 4
set key font "roman, 12"  samplen 1  

set grid

plot \
     'zgoubi.OPTICS.out' u ($13/100.):15  w l lc 1  tit "Y orbit"  ,\
     'zgoubi.OPTICS.out' u ($13/100.):17  w l lc 1  tit "Z orbit"

set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_OPTICS_YZorbits.eps" 
 replot 
 set terminal X11 
 unset output 

pause 3


