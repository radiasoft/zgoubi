
      set title "Tune vs turn number"  font "roman,16" 

      set xlabel "Q_x, Q_y" font "roman,12" 
      set ylabel "x-amplitude" font "roman,12" 
      set y2label "y-amplitude" font "roman,12" 

      set xtics mirror font "roman,11" 
      set ytics nomirror font "roman,11" 
      set y2tics nomirror font "roman,11"  

       set key maxcol 1
       set key c t

      set logscale y
      set logscale y2

       plot  \
       "tunesFromFai_spctra.Out"   u ($3):($4) axes x1y2 w lp ps .1 linecolor rgb "blue"  tit "Q_y" ,\
       "tunesFromFai_spctra.Out"   u ($1):($2) axes x1y1 w lp ps .1 linecolor rgb "red" tit "Q_x"  


      set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_spectrumFromFai_xy.eps"  
       replot  
       set terminal X11  
       unset output  

      pause 1

      set xlabel "Q_l [{/Symbol \264}10^{-3}]" font "roman,12" 
      set ylabel "l-amplitude" font "roman,14" 
      unset y2label 

      unset logscale y
      unset logscale y2

      set xtics font "roman,12" 
      set ytics mirror font "roman,12" 
      unset y2tics 

      plot  \
       "tunesFromFai_spctra.Out"   u ($5 *1e3):($6)  w lp ps .3 linecolor rgb "blue"  tit "Q_l" 

      set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_spectrumFromFai_l.eps"  
       replot  
       set terminal X11  
       unset output  

      pause 8


exit
