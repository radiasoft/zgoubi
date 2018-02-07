      set title "Tunes versus turn #"  font "roman,16" 

      set xlabel "Turn # (*1000)" font "roman,14" 
      set ylabel "Qx, Qy" font "roman,14" 
#      set y2label "mQx - nQy" font "roman,14" 

      set xtics     mirror font "roman,11" 
      set ytics   mirror font "roman,11" 
#      set y2tics  nomirror font "roman,11" 

      set key font "roman,15" center center
      set grid

 G = 1.1596e-3
 m = 1 
 n = -2

      plot \
       "tunesFromFai_iterate.out"   u ($15 /1000):($3) axis x1y1 w lp  pt 4 ps .2 linecolor rgb "red" tit "Qx" ,\
       "tunesFromFai_iterate.out"   u ($15 /1000):($4) axis x1y1 w lp  pt 4 ps .2 linecolor rgb "blue" tit "Qy"

       set terminal postscript eps blacktext color  enh  "Times-Roman" 14  
       set output "gnuplot_QxQy.vs.turn.eps"  
       replot  
       set terminal X11  
       unset output  


      pause 1


      set title "Emittance (or invariant) versus turn #"  font "roman,16" 

      set xlabel "Turn # (*1000)" font "roman,14"  
      set ylabel  "{/Symbol e}_x/{/Symbol p} [m]" font "roman,14"  
      set y2label  "{/Symbol e}_y/{/Symbol p} [m]" font "roman,14"  

      set xtics           font "roman,11"            
      set ytics  nomirror font "roman,11"  
      set y2tics nomirror font "roman,11" 

      set key font "roman,15" top left
      set grid

 set logscale y ; set format y "%.0s*10^{%T}"
 set logscale y2 ; set format y2 "%.0s*10^{%T}"

      plot [][]\
       "tunesFromFai_iterate.out"   u ($15 /1000):($7) axes x1y1 w lp  pt 4 ps .2 linecolor rgb "red" tit "{/Symbol e}_x" ,\
       "tunesFromFai_iterate.out"   u ($15 /1000):($8) axes x1y2 w lp  pt 4 ps .2 linecolor rgb "blue"  tit "{/Symbol e}_y" 
#,\
#       "tunesFromFai_iterate.out"   u ($15 /1000):($7/m + $8/n) axis x1y2 w lp pt 5 ps .8  tit "ex/m-ey/n" 

       set terminal postscript eps blacktext color  enh  "Times-Roman" 14  
       set output "gnuplot_exey.vs.turn.eps"  
       replot  
       set terminal X11  
       unset output  


      pause 8

 exit
