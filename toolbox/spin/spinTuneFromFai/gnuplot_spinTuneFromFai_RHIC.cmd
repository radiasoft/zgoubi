      #Tunes versus turn number !  
      set xlabel "TURN REGION" font "roman,28" 
      set ylabel "Q_s" font "roman,28" 

      set title "TUNE / CYCLE"  font "roman,28" 
      set xtics font "roman,18" 
      set ytics nomirror font "roman,18" 

#      set yrange [.48:.52]

      plot \
       "spinTuneFromFai_iterate.Out"   u ($5):($1) w p ps .8 tit "Q_s" 

      set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_Qspin.vs.turn.eps"  
       replot  
       set terminal X11  
       unset output  

      pause 80


 exit
