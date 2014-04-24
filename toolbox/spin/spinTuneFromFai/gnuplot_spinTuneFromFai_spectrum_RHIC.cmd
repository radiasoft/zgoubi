      #Tunes versus turn number !  
      set xlabel "Frac. Q_s" font "roman,28" 
      set ylabel "Amplitude" font "roman,28" 

      set title "Spin tune spectra (superimposed)"  font "roman,28" 
      set xtics font "roman,18" 
      set ytics nomirror font "roman,18" 

#      set xrange [0.:1.]

      plot \
       "spinTuneFromFai_spectrum_iterate.Out"   u ($1):($2) w l  notit 
#       "spinTuneFromFai_spectrum.Out"   u ($1):($2) w l  notit 

      pause 8

      set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_QspinSpectrum.eps"  
       replot  
       set terminal X11  
       unset output  

 exit
