      #Tunes versus turn number !  

      set title "Emittances, from tunesFromFai"  font "roman,16" 

      set xlabel "E [MeV]" font "roman,16" 
      set ylabel  "{/Symbol e}_{x,y}/{/Symbol p}, norm." font "roman,16" 
#      set y2label "{/Symbol e}_y/{/Symbol p}, norm." font "roman,16" 

      set xtics font "roman,12" 
      set ytics mirror font "roman,12" 
#      set y2tics nomirror font "roman,12"  

#      set yrange [.6848:.6852]
#      set y2range [.6748:.6752]

      set key  t l
     set key maxrow 1

       set logscale y
       set logscale y2
       set format y "10^{%L}"
       set format y2 "10^{%L}"

am = 0.511

      plot [40:165] \
       "tunesFromFai.out_H"   u ($11):($7 *$11/am) axes x1y2 w l  tit "{/Symbol e}_x H" ,\
       "tunesFromFai.out_V"   u ($11):($8 *$11/am) axes x1y1 w l smooth cspline tit "{/Symbol e}_y V" 

#       "tunesFromFai.out_Hz"   u ($11):($7 *$11/am) axes x1y1 w l  tit "{/Symbol e}_x Hz" ,\

      set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_tunesFromFai_emittances.eps"  
       replot  
       set terminal X11  
       unset output  

      pause 2

 exit
