      #Tunes versus turn number !  
      set xlabel "Turn #" font "roman,14" 
      set ylabel "Qx, Qy" font "roman,14" 
#      set y2label "Dist. to resonance" font "roman,14" 

      set title "Tunes versus turn #"  font "roman,14" 

      set xtics font "roman,16" 
      set ytics nomirror font "roman,14" 
      set y2tics nomirror font "roman,14"  

      set key center top
      set grid

      plot \
       "tunesFromFai.out"   u ($15):($3) axis x1y1 w lp  ps .8  tit "Qx" ,\
       "tunesFromFai"   u ($15):($4) axis x1y1 w lp  ps .8  tit "Qy" 

      set terminal postscript eps blacktext color  enh size 8.3cm,4.4cm "Times-Roman" 14  
       set output "gnuplot_Tunes.vs.Ekin.eps"  
       replot  
       set terminal X11  
       unset output  


      pause 8

      plot \
       "tunesFromFai_iterate.out"   u ($15):($7) axis x1y1 w lp pt 4 ps 1.8  tit "x emittance" ,\
       "tunesFromFai_iterate.out"   u ($15):($7) axis x1y1 w lp pt 5 ps 1.8  tit "y emittance" 

      set terminal postscript eps blacktext color  enh size 8.3cm,4.4cm "Times-Roman" 14  
       set output "gnuplot_exey.vs.turn.eps"  
       replot  
       set terminal X11  
       unset output  


      pause 888

 exit
