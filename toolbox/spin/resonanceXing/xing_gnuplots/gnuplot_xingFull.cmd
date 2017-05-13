 set xlabel "Turn"       font "roman,20"
 set ylabel "Sy "               font "roman,20"
 set title "RESONANCE CROSSING" font "roman,20"

 set xtics font "roman,14"
 set ytics font "roman,14"

#plot full xing.  17-19 : SX,-Y,-Z,|S|
#plot 'fromBFai2Fai.out' us ($24):($22)  w p pt 3 ps .1   # SZ vs. turn#
#plot 'fromBFai2Fai.out' us ($24):($22)  w p pt 3 ps .1  title 'Zgoubi'    # SZ vs. kin-energy 
 plot \
  'fromBFai2Fai.out' us ($24):($38<180 ?   $22 : 1/0)  w l lt 3 lw .5  title 'Zgoubi' 


 set terminal postscript eps blacktext color
 set output "gnuplot_xingFull.eps"
 replot
 set terminal X11
 unset output

pause .1

exit


