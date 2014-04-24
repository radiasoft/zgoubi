 set xlabel "x (m)"       font "roman,20"
 set ylabel "x' (rad) "               font "roman,20"
 set title "X-X' PHASE SPACE" font "roman,18"

 set xtics font "roman,14"
 set ytics font "roman,14"

 plot 'fromBFai2Fai.out' us ($10/100.):($11/1000.)  w p pt 3 ps .4  title 'Zgoubi'   

#pause .1 

 set terminal postscript eps blacktext color
 set output "gnuplot_xing-xxp.eps"
 replot
 set terminal X11
 unset output

exit


