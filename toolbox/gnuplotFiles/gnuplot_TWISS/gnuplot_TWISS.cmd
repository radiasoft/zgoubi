
 set title "Optical functions, from zgoubi.TWISS.out"    font "roman,16"   # offset 0,+.7    

 set xlabel "s [m]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "{/Symbol b}_x, {/Symbol b}_y [m]"             font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "{/Symbol h}_x, {/Symbol h}_y"          font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set x2tics  font "roman,12" mirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics  font "roman,12" nomirror      #offset 0,-.6

#set key above maxrows 1 
#set key tm c maxrows 1 width 4
set key t c maxrows 1 width 4
set key font "roman, 12"  samplen 1  

set grid

set samples 10000
set xrange  []
set x2range []
set yrange []
set y2range []

set xrange  []

#set xrange  [490:510]

 plot \
      "zgoubi.TWISS.out" u ($13):($2) axes x1y1 w l lt 1 lc rgb "red" tit "{/Symbol b}_x"  ,\
      "zgoubi.TWISS.out" u ($13):($4) axes x1y1 w l lt 1 lc rgb "blue" tit "{/Symbol b}_y"  ,\
      "zgoubi.TWISS.out" u ($13):($7) axes x1y2 w l lt 1 lc "magenta" lw 2 tit "{/Symbol h}_x"  ,\
      "zgoubi.TWISS.out" u ($13):($9) axes x1y2 w l lt 1 lc rgb "cyan" lw 2 tit "{/Symbol h}_y"  

 set samples 10000
 set terminal postscript eps blacktext color enh size 8.3cm,5cm "Times-Roman" 12
 set output "gnuplot_TWISS_btxy.eps"
 replot
 set terminal X11
 unset output

pause 2

 set title "Orbit, from zgoubi.TWISS.out"    font "roman,16"   # offset 0,+.7    

 set ylabel "x, y [m]"             font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "x, y [m]"             font "roman,13"   #offset -0,-1 rotate by -20 

 unset x2tics 
 set xtics  font "roman,12" mirror
 unset y2tics
 set ytics  font "roman,12" mirror      #offset 0,-.6

 plot \
      "zgoubi.TWISS.out" u ($13):($15) axes x1y1 w l lt 1 lc rgb "red" tit "x"  ,\
      "zgoubi.TWISS.out" u ($13):($17) axes x1y1 w l lt 1 lc rgb "blue" tit "y"  

 set samples 10000
 set terminal postscript eps blacktext color enh size 9.3cm,6cm "Times-Roman" 12
 set output "gnuplot_TWISS_xy.eps"
 replot
 set terminal X11
 unset output

      pause 2

 exit
