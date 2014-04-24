      #Tunes versus turn number !  
      set xlabel "TURN REGION" font "roman,28" 
      set ylabel "Qx" font "roman,28" 
      set y2label "Qy" font "roman,28" 

      set title "TUNE / CYCLE"  font "roman,28" 
      set xtics font "roman,18" 
      set ytics nomirror font "roman,18" 
      set y2tics nomirror font "roman,18"  

#      set yrange [.6848:.6852]
#      set y2range [.6748:.6752]

      plot \
       "tunesFromFai_iterate.out"   u ($15):($3>.55 && $3<.85 ? $3 : 1/0) axis x1y1 w p pt 1 ps 1.1  tit "Qx" ,\
       "tunesFromFai_iterate.out"   u ($15):($5>.55 && $5<.85 ? $5 : 1/0) axis x1y1 w p pt 1 ps 1.1  tit "Qx" ,\
       "tunesFromFai_iterate.out"   u ($15):($4>.55 && $4<.85 ? $4 : 1/0) axis x1y2 w p pt 2 ps 1.1  tit "Qy" ,\
       "tunesFromFai_iterate.out"   u ($15):($6>.55 && $6<.85 ? $6 : 1/0) axis x1y2 w p pt 2 ps 1.1  tit "Qy"


      set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_QxQy.vs.turn.eps"  
       replot  
       set terminal X11  
       unset output  

      pause 88

 exit
