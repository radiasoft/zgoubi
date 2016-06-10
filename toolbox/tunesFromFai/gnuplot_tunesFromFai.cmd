      #Tunes versus turn number !  

      set title "TUNE / CYCLE"  font "roman,16" 

      set xlabel "TURN REGION" font "roman,16" 
      set ylabel "Qx" font "roman,16" 
      set y2label "Qy" font "roman,16" 

      set xtics font "roman,12" 
      set ytics nomirror font "roman,12" 
      set y2tics nomirror font "roman,12"  

#      set yrange [.6848:.6852]
#      set y2range [.6748:.6752]

      plot \
       "tunesFromFai.out"   u ($11):($3>.51 ? $5 : 1/0) axis x1y1 w p pt 1 ps 1.1  tit "Qx" ,\
       "tunesFromFai.out"   u ($11):($5>.51 ? $3 : 1/0) axis x1y1 w p pt 1 ps 1.1  tit "Qx" ,\
       "tunesFromFai.out"   u ($11):($4>.51 ? $6 : 1/0) axis x1y2 w p pt 2 ps 1.1  tit "Qy" ,\
       "tunesFromFai.out"   u ($11):($6>.51 ? $4 : 1/0) axis x1y2 w p pt 2 ps 1.1  tit "Qy"

      set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_tunesFromFai.eps"  
       replot  
       set terminal X11  
       unset output  

      pause 88

 exit
