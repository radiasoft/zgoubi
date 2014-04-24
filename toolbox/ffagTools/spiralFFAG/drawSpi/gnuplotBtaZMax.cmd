 
set xlabel 'Spiral angle xi'
set ylabel 'Beta_Z max'
set title 'Maximum beta_Z amplitude '

### plot [0:80] [0:15] "scanKXi.out" using ($3):($7) with lines  
 plot [] [0:15] "scanKXi.out" using ($3):($9<998 ? $7 : 1/0)  with p pt 6 ps 1.4 title "Parameter is K", "scanKXi.out" using ($3):($9>998 ? $7 : 1/0) with points pointtype 7 pointsize 2.5 title "Sample tuning"

##set multiplot
## plot [0:70] [0:15] "scanKXi.out" using ($3):($9<99 ? $7 : 1/0)  with lp pt 6 ps 1.4 
## plot [0:70] [0:15]  "scanKXi.out" using ($3):($9>998 ? $7 : 1/0) w p pointtype 7 pointsize 2
##unset multiplot

pause 2

 set terminal postscript eps blacktext color
 set output "gnuplotBtaZMax.eps"
 replot

set terminal X11
set output

