
 set title "Orbits, from zgoubi.OPTICS.out"    font "roman,16"   # offset 0,+.7    

 set xlabel "s [m]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "x, y [m]"             font "roman,13"   #offset -0,-1 rotate by -20 
# set y2label "{/Symbol h}_x, {/Symbol h}_y"          font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12"   mirror
# set x2tics  font "roman,12" mirror
 set ytics  font "roman,12"   mirror      #offset 0,-.6
# set y2tics  font "roman,12" nomirror      #offset 0,-.6

set key t l maxrows 1 width 4
set key font "roman, 12"  samplen 1  

Circ = 3.88294E+05/100.

set grid

plot \
     'zgoubi.OPTICS.out' u ($13/100. + ($30 -1)*Circ):15  w l lc 1  tit "x orbit"  ,\
     'zgoubi.OPTICS.out' u ($13/100. + ($30 -1)*Circ):17  w l lc 3  tit "y orbit"

set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_OPTICS_xy.eps" 
 replot 
 set terminal X11 
 unset output 

pause 3

 set title "{/Symbol b}_{x,y}, {/Symbol h}_{x,y}, from zgoubi.OPTICS.out"    font "roman,16"   # offset 0,+.7    

 set xlabel "s [m]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "{/Symbol b}_x, {/Symbol b}_y  [m]"             font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "{/Symbol h}_x, {/Symbol h}_y [m]"          font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12"   
# set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror  
 set y2tics  font "roman,12" nomirror      #offset 0,-.6

set key t c maxrows 1 width 4
set key font "roman, 12"  samplen 1  

set grid

plot \
     'zgoubi.OPTICS.out' u ($13/100. + ($30 -1)*Circ):2  w l lc 1  tit "{/Symbol b}_x"  ,\
     'zgoubi.OPTICS.out' u ($13/100. + ($30 -1)*Circ):4  w l lc 3  tit "{/Symbol b}_y" ,\
     'zgoubi.OPTICS.out' u ($13/100. + ($30 -1)*Circ):7  axes x1y2 w l lc 1  tit "{/Symbol h}_x"  ,\
     'zgoubi.OPTICS.out' u ($13/100. + ($30 -1)*Circ):9  axes x1y2 w l lc 3  tit "{/Symbol h}_y"

set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_OPTICS_bxyDxy.eps" 
 replot 
 set terminal X11 
 unset output 

pause 3


 set title "{/Symbol b}_{x,y}, {/Symbol h}_{x,y}, from zgoubi.OPTICS.out"    font "roman,16"   # offset 0,+.7    

 set xlabel "s [m]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "{/Symbol b}_x, {/Symbol b}_y  [m]"             font "roman,13"   #offset -0,-1 rotate by -20 
 unset y2label #"{/Symbol h}_x, {/Symbol h}_y [m]"          font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12"   
# set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror  
 unset y2tics # font "roman,12" nomirror      #offset 0,-.6

set key t c maxrows 1 width 4
set key font "roman, 12"  samplen 1  

set grid

plot \
     'zgoubi.OPTICS.out' u ($13/100. + ($30 -1)*Circ):2  w l lc 1  tit "{/Symbol b}_x"  ,\
     'zgoubi.OPTICS.out' u ($13/100. + ($30 -1)*Circ):4  w l lc 3  tit "{/Symbol b}_y" 

set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_OPTICS_bxy.eps" 
 replot 
 set terminal X11 
 unset output 

pause 3


