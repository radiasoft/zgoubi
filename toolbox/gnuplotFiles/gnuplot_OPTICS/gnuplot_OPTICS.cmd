
plot [1137:1140.5] 'zgoubi.OPTICS.out' u ($13/100.):2  w l lc 1 , './madFiles/madzg.in' u 12:14 w l lc 3

set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_OPTICS.eps" 
 replot 
 set terminal X11 
 unset output 

pause 3


