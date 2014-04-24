 set xlabel "phase (rad)"       font "roman,20"
 set ylabel "dp/p "               font "roman,20"
 set title "l-dp PHASE SPACE" font "roman,18"

 set xtics font "roman,14"
 set ytics font "roman,14"

 plot 'fromBFai2Fai.out' us ($34):($35)  w p pt 3 ps .4  title 'Zgoubi'   

#pause .1

 set terminal postscript eps blacktext color
 set output "gnuplot_xing-ldp.eps"
 replot
 set terminal X11
 unset output

exit




