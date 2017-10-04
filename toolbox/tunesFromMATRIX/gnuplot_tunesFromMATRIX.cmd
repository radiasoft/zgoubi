

 set title "Tunes"    font "roman,16"   # offset 0,+.7    

 set xlabel "E [MeV]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "Q_x, Q_y"             font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "x_{orbit} (mm)"          font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics  font "roman,12" nomirror      #offset 0,-.6

set key bot left  font "roman, 12"  samplen 1  

set grid

 plot \
      "tunesFromMatrix.out" u ($12):($3) axes x1y1 w lp ps .2 lt 20 tit "Q_x"  ,\
      "tunesFromMatrix.out" u ($12):($4) axes x1y1 w lp ps .2 lt 30 tit "Q_y"  ,\
      "tunesFromMatrix.out" u ($12):($1*1e3) axes x1y2 w lp ps .2 lt 40 tit "Orbit"  

set samples 100000
 set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12
 set output "gnuplot_tunesFromMATRIX_QxQy.eps"
 replot
 set terminal X11
 unset output

      pause 1

 exit