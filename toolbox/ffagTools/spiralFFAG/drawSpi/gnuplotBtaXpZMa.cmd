 
set xlabel 'Field index K'
set ylabel 'spiral angle (deg)'
set title 'Stable K/xi diagram'

 plot [0:7] [0:70] "scanKXi.out" using 2:3 with points pointtype 4 pointsize 1.5 

pause 2

 set terminal postscript eps blacktext color
 set output "gnuplotKxi.eps"
 replot
set terminal X11
set output

#set multiplot        #  will cause all following plots to be plotted simultaneously
# plot "scanKXi.out" using 3:6 with lines
#pause 1
# plot "scanKXi.out" using 3:7 with lines
#pause 100
#unset multiplot
