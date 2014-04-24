
set tit 'beta_x (m) versus s (m)'
set xlab 's (m)'
set ylab 'beta_x (m)'

plot [0:808] [] 'betaFromMatrix.out' us ($1/100.):3 w l lw 2  tit 'zgoubi' 

pause 2

 set terminal postscript eps blacktext color
 set output "gnuplot_betx.eps"
 replot
 set terminal X11
 set output



set tit 'beta_y (m) versus s (m)'
set xlab 's (m)'
set ylab 'beta_y (m)'

plot [0:808] [] 'betaFromMatrix.out' us ($1/100.):7 w l lw 2  tit 'zgoubi'

pause 2

 set terminal postscript eps blacktext color
 set output "gnuplot_bety.eps"
 replot
 set terminal X11
 set output


set tit 'D_x (m) versus s (m)'
set xlab 's (m)'
set ylab 'D_x (m)'

plot [0:808] [] 'betaFromMatrix.out' us ($1/100.):4 w l lw 2  tit 'zgoubi'

pause 2

 set terminal postscript eps blacktext color
 set output "gnuplot_Dx.eps"
 replot
 set terminal X11
 set output



