
 set title "Tunes and chromas from Twiss scan"    font "roman,16"   # offset 0,+.7    

 set xlabel "s [m]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "{/Symbol b}_x, {/Symbol b}_y [m]"             font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "{/Symbol h}_x, {/Symbol h}_y"          font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set x2tics  font "roman,12" mirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics  font "roman,12" nomirror      #offset 0,-.6

#set key above maxrows 1 
set key tm c maxrows 1 width 4
set key font "roman, 12"  samplen 1  

set grid

set samples 10000
set xrange  []
set x2range []
set yrange []
set y2range []

 plot \
      "zgoubi.TWISS.out" u ($13/100.):($2) axes x1y1 w l lt 1 lc 1 tit "{/Symbol b}_x"  ,\
      "zgoubi.TWISS.out" u ($13/100.):($4) axes x1y1 w l lt 1 lc 3 tit "{/Symbol b}_y"  ,\
      "zgoubi.TWISS.out" u ($13/100.):($7) axes x1y2 w l lt 1 lc 2 lw 2 tit "{/Symbol h}_x"  ,\
      "zgoubi.TWISS.out" u ($13/100.):($9) axes x1y2 w l lt 1 lc 4 lw 2 tit "{/Symbol h}_y"  

 set samples 10000
 set terminal postscript eps blacktext color enh size 8.3cm,6cm "Times-Roman" 12
 set output "gnuplot_TWISS.eps"
 replot
 set terminal X11
 unset output

      pause 1

 exit
