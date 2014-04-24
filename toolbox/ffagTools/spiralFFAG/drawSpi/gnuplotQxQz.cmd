 
set xlabel 'Qx'
set ylabel 'Qz'
set title 'Tune diagram' 

set grid

# plot [0.:.5] [0.:.5] "scanKXi.out" using 4:5 with points pointtype 4 pointsize 1.5 
# plot  [0.:5] [0.:5] "scanKXi.out" using 4:($9<999 ? $5 : 1/0)  with points pointtype 4 pointsize .8  title "Stable motion", "scanKXi.out" using 4:($9>998 ? $5 : 1/0)  with points pointtype 5 pointsize 2  title "Sample tuning" 
plot  [0.:10] [0.:10] "scanKXi.out" using 4:($9<999 ? $5 : 1/0)  with points pointtype 4 pointsize .8  title "Stable motion", "scanKXi.out" using 4:($9>998 ? $5 : 1/0)  with points pointtype 5 pointsize 2  title "Sample tuning"

pause 2

 set terminal postscript eps blacktext color
 set output "gnuplotQxQz.eps"
 replot
set terminal X11
set output

#set multiplot        #  will cause all following plots to be plotted simultaneously
# plot "scanKXi.out" using 4:5 with dots
#pause 100
#unset multiplot

unset grid 
