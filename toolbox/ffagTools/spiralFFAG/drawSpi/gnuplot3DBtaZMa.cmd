 
set xlabel 'Field index K'
set ylabel 'Spiral angle xi'
set title 'Maximum beta_Z amplitude '

set dgrid3d
set hidden3d 
set pm3d

 splot [0:7] [0:70]  [0:15] "scanKXi.out" using 2:3:7 with lines 

pause 2

 set terminal postscript eps blacktext color
 set output "gnuplot3DBtaZMa.eps"
 replot
set terminal X11
set output

unset hidden3d
unset pm3d
unset dgrid3d
