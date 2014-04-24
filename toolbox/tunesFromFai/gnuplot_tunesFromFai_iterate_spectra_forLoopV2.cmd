      #Tunes versus turn number !  
      set xlabel "TURN REGION" font "roman,28" 
      set ylabel "Qx" font "roman,28" 
      set y2label "Qy" font "roman,28" 

      set title "TUNE / CYCLE"  font "roman,28" 
      set xtics font "roman,18" 
      set ytics nomirror font "roman,18" 
      set y2tics nomirror font "roman,18"  


#      plot [4.7:5.3][] for [var = 1 : 300 : 200] \
#       "tunesFromFai_iterate_spctra.Out"   u ($1):($7==var ? $2 : 1/0)  w p ps .3  tit "Qx"  

#      plot [][] for [var = 1001 : 1601 : 200] \
#       "tunesFromFai_iterate_spctra.Out"   u ($3):($7==var && $1>5.5 ? $4 : 1/0)  w p ps .3  tit "Qy"  

      plot [][] for [var = 2401 : 2801 : 200] \
       "tunesFromFai_iterate_spctra.Out"   u ($1):($7==var && $1>4.5 && $1<5 ? $2 : 1/0)  w p ps .3  tit "Qx"  

#,\
#       "tunesFromFai_iterate_spctra.Out"   u ($1):($7==1 ? $4 : 1/0)  w p ps .3  tit "Qy" 

      pause 888

      set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_spectrum.eps"  
       replot  
       set terminal X11  
       unset output  

 exit
