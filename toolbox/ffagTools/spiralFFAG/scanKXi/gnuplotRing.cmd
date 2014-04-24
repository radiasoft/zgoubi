 
set title 'Spiral ring '
set size ratio -1

 plot "drawSpi.out"  with lines  

pause 2

 set terminal postscript eps blacktext color
 set output "gnuplotRing.eps"
 replot
set terminal X11
set output


