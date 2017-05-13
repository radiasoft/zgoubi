 set xlabel "Turn"       font "roman,20"
 set ylabel "Sy "               font "roman,20"
 set title "RESONANCE CROSSING, FIRST 30 TURNS" font "roman,18"

 set xtics font "roman,14"
 set ytics font "roman,14"

#plot full xing.  17-19 : SX,-Y,-Z,|S|
#plot 'fromBFai2Fai.out' us ($38):($22)  w p pt 3 ps .1   # SZ vs. turn#
# plot 'fromBFai2Fai.out' us ($24):($38<200 ?   $22 : 1/0)  w l lt 3 lw .5  title 'Zgoubi' 

  plot 'fromBFai2Fai.out' us ($38):($38<50 ?   $22 : 1/0)  w l lt 3 lw .5  title 'Zgoubi' 

#pause .1 

 set terminal postscript eps blacktext color
 set output "gnuplot_xingInit.eps"
 replot
 set terminal X11
 unset output

exit


