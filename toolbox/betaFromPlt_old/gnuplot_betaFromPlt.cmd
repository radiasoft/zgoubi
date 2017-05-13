
 set title "Beta and dispersion functions \n - from betaFromPlt -"    font "roman,16"   # offset 0,+.7    

 set xlabel "s (m)"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "{/Symbol b}_x, {/Symbol b}_y (m)"             font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "D_x, D_y (m)"          font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics  font "roman,12" nomirror      #offset 0,-.6

set key bot left  font "roman, 12"  samplen 1  
set key maxrow 1
set grid

 plot  \
      "betaFromPlt.out" u ($13):($2)  axes x1y1 w l lt 1 lc 1 tit "{/Symbol b}_x"  ,\
      "betaFromPlt.out" u ($13):($4)  axes x1y1 w l lt 1 lc 3 tit "  {/Symbol b}_y"  ,\
      "betaFromPlt.out" u ($13):($7)  axes x1y2 w l lt 2 lc 1 tit "  D_x"  ,\
      "betaFromPlt.out" u ($13):($9)  axes x1y2 w l lt 2 lc 3 tit "  D_y"  

#set samples 100000
 set terminal postscript eps blacktext color enh size 8.3cm,6cm "Times-Roman" 12
 set output "gnuplot_betaFromPlt.eps"
 replot
 set terminal X11
 unset output

      pause 1

 exit
