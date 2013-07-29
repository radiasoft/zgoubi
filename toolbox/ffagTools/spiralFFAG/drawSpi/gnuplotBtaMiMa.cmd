 
set xlabel 'beta_X (m)'
set ylabel 'bate_Z (m)'
set title '  '

### plot [1:8.] [2:10.] "scanKXi.out" using 6:7 with linespoints pointtype 6 pointsize 1.4 
plot  "scanKXi.out" using 6:($9<999 ? $7 : 1/0)  with linespoints pointtype 6 pointsize 1.4 title "K, xi scan","scanKXi.out" using 6:($9>998 ? $7 : 1/0)  with points pointtype 7 pointsize 2.5 title "Sample tuning"

pause 2

 set terminal postscript eps blacktext color
 set output "gnuplotBtaMiMa.eps"
 replot
set terminal X11
set output


