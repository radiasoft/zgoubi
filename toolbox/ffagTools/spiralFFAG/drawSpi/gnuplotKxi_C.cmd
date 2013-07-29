 
set xlabel 'Field index K'
set ylabel 'spiral angle (deg)'
set title 'Stable K/xi diagram'

set grid 


# carbon scales
 plot [0:30] [0:90] "scanKXi.out" using 2:($9<999 ? $3 : 1/0) with points pointtype 4 pointsize 1.5 title "Stable motion", "scanKXi.out" using 2:($9>998 ? $3 : 1/0) with points pointtype 5 pointsize 2.5 title "Sample tuning"


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

unset grid 
