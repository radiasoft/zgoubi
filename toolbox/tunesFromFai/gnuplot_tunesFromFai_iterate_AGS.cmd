      set title "Tunes versus turn #"  font "roman,16" 

      set xlabel "Turn # (*1000)" font "roman,14" 
      set ylabel "Qx, Qy" font "roman,14" 
      set y2label "mQx - nQy" font "roman,14" 

      set xtics     mirror font "roman,11" 
      set ytics   nomirror font "roman,11" 
      set y2tics  nomirror font "roman,11" 

      set key font "roman,15" center center
      set grid

 G = 1.7928474/938.27203
 m = 1 
 n = -2

      plot \
       "tunesFromFai_iterate.out"   u ($15 /1000):($3) axis x1y1 w lp  ps .8  tit "Qx" ,\
       "tunesFromFai_iterate.out"   u ($15 /1000):($4) axis x1y1 w lp  ps .8  tit "Qy" ,\
       "tunesFromFai_iterate.out"   u ($15 /1000):(m*$3+n*$4) axis x1y2 w lp  ps .6  tit "mQx+nQy"

       set terminal postscript eps blacktext color  enh size 8.3cm,4.4cm "Times-Roman" 14  
       set output "gnuplot_QxQy.vs.turn.eps"  
       replot  
       set terminal X11  
       unset output  


      pause 1


      set title "Emittances versus turn #"  font "roman,16" 

      set xlabel "Turn # (*1000)" font "roman,14"  
      set ylabel  "emittance (pi.m.rd)" font "roman,14"  
      set y2label "  ex/m + ey/n" font "roman,14" 

      set xtics           font "roman,11"            
      set ytics  nomirror font "roman,11"  
      set y2tics nomirror font "roman,11" 

      set key font "roman,15" top left
      set grid

      plot [][]\
       "tunesFromFai_iterate.out"   u ($15 /1000):($7) axis x1y1 w lp pt 4 ps .8  tit "x emit." ,\
       "tunesFromFai_iterate.out"   u ($15 /1000):($8) axis x1y1 w lp pt 5 ps .8  tit "y emit." ,\
       "tunesFromFai_iterate.out"   u ($15 /1000):($7/m + $8/n) axis x1y2 w lp pt 5 ps .8  tit "ex/m-ey/n" 

       set terminal postscript eps blacktext color  enh size 8.3cm,4.4cm "Times-Roman" 14  
       set output "gnuplot_exey.vs.turn.eps"  
       replot  
       set terminal X11  
       unset output  


      pause 888

 exit
