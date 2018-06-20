
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

Circ = 0. 

set grid

#plot [:7300] [-.04:.04] \

plot [:] [:] \
     'zgoubi.OPTICS.out'  u ($13 + ($30 -1)*Circ):15  w l lc 1  tit "x orbit"  ,\
     'zgoubi.OPTICS.out'  u ($13 + ($30 -1)*Circ):17  w l lc 3  tit "y orbit"

set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_OPTICS_xy.eps" 
 replot 
 set terminal X11 
 unset output 

pause .5

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
 
# "every 3" in this plot assumes there are 3 different label classes stored in zgoubi.OPTICS.out
plot [:] \
     'zgoubi.OPTICS.out' every 1  u ($13 + ($30 -1)*Circ):7  axes x1y2 w l lc 10 lt 2  tit "{/Symbol h}_x"  ,\
     'zgoubi.OPTICS.out' every 1  u ($13 + ($30 -1)*Circ):9  axes x1y2 w l lc 30 lt 2 tit "{/Symbol h}_y" ,\
     'zgoubi.OPTICS.out' every 1  u ($13 + ($30 -1)*Circ):($2)  w l lc 1  lt 1 tit "{/Symbol b}_x"  ,\
     'zgoubi.OPTICS.out' every 1  u ($13 + ($30 -1)*Circ):($4)  w l lc 3  lt 1 tit "{/Symbol b}_y" 

set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_OPTICS_bxyDxy.eps" 
 replot 
 set terminal X11 
 unset output 

pause .5

 set title "H & V determinants and momentum \n from zgoubi.OPTICS.out"    font "roman,16"   # offset 0,+.7    

 set xlabel "s [m]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "Y, Z determinants"             font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "p [MeV/c]"          font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12"   
# set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror  
 set y2tics  font "roman,12" nomirror      #offset 0,-.6

set key t c maxrows 1 width 4
set key font "roman, 12"  samplen 1  
set key maxrow 1

set grid

p0 = 123.

#plot [2000:2250] \

plot [:] \
     'zgoubi.OPTICS.out' u ($13/100. + ($30 -1)*Circ):(1.-($33*$36-$34*$35))  w l lc rgb "red"  tit '1-Det_Y'  ,\
     'zgoubi.OPTICS.out' u ($13/100. + ($30 -1)*Circ):(1.-($37*$40-$38*$39))  w l lc rgb "blue"  tit '1-Det_Z' ,\
     'zgoubi.OPTICS.out' u ($13/100. + ($30 -1)*Circ):(($31+$32) *p0) axes x1y2 w l lc rgb "cyan"  tit 'Momentum'

set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_OPTICS_DetYZ.eps" 
 replot 
 set terminal X11 
 unset output 

pause .5

exit 
