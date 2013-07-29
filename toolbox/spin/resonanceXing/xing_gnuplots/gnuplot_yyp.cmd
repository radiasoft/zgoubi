 set xlabel "y (m)"       font "roman,20"
 set ylabel "y' (rad) "               font "roman,20"
 set title "Y-Y' PHASE SPACE" font "roman,18"

 set xtics font "roman,14"
 set ytics font "roman,14"

 plot 'fromBFai2Fai.out' us ($12/100.):($13/1000.)  w p pt 3 ps .4  title 'Zgoubi'  

 set terminal postscript eps blacktext color
 set output "gnuplot_xing-yyp.eps"
 replot
 set terminal X11
 unset output

#pause .1 

exit


