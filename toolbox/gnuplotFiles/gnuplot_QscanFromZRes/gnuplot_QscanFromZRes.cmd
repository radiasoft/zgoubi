#A: line increment
#B: data block increment
#C: The first line
#D: The first data block
#E: The last line
#F: The last data block


system "grep 'NU_Y = ' zgoubi.res | cat > gnuplot_temp"

 set title "Optical functions, from zgoubi.TWISS.out"    font "roman,16"   # offset 0,+.7    

 set xlabel "iteration number"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "Q1, Q2"             font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6

#set key above maxrows 1 
#set key tm c maxrows 1 width 4
set key t c maxrows 1 width 4
set key font "roman, 12"  samplen 1  

set grid

# skip line # 0

 plot \
      "gnuplot_temp" every ::1 u :($3) w p pt 4 ps .4 lc rgb "red" tit "Q_1"  ,\
      "gnuplot_temp" every ::1 u :($6) w p pt 5 ps .4 lc rgb "blue" tit "Q_2"  

 set samples 10000
 set terminal postscript eps blacktext color enh size 8.3cm,5cm "Times-Roman" 12
 set output "gnuplot_Qscan.eps"
 replot
 set terminal X11
 unset output

pause 1
exit
