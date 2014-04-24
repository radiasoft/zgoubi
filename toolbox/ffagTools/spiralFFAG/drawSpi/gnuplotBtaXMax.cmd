 
set xlabel 'Index K'
set ylabel 'Beta_X max'
set title 'Maximum beta_X amplitude '

# plot [0:10] [0:15] "scanKXi.out" using 2:6 with lines  
 plot  [] [0:15]  "scanKXi.out" using 2:($9<999 ? $6 : 1/0) with points pointtype 8 pointsize 1.6 title "Parameter is xi",  "scanKXi.out" using 2:($9>998 ? $6 : 1/0) with p pointtype 9 pointsize 2.5 title "Sample tuning"

pause 2

 set terminal postscript eps blacktext color
 set output "gnuplotBtaXMax.eps"
 replot

set terminal X11
set output

