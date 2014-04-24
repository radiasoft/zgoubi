 set xtics font "roman,18"
 set ytics font "roman,18"
 
set grid
 
set tit  "beta_x,y versus s" font "roman,25"
set xlab "s (m)"             font "roman,25"
set ylab "beta_x, -beta_y (m)"  font "roman,25"
plot \
 "betaFromMatrix.out" us ($1/100.):($3) w l lw 2  tit "beta_x" , \
 "betaFromMatrix.out" us ($1/100.):($7 * (-1)) w l lw 2  tit "-beta_y" 
 
 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_betxy.eps"
 replot
 set terminal X11
 unset output
 
 
set tit  "D_x, D_y versus s" font "roman,25"
set xlab "s (m)"             font "roman,25"
set ylab "D_x, D_y (m)"      font "roman,25"
plot \
 "betaFromMatrix.out" us ($1/100.):($4) w l lw 2  tit "D_x" , \
 "betaFromMatrix.out" us ($1/100.):($8) w l lw 2  tit "D_y" 
 
 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_Dxy.eps"
 replot
 set terminal X11
 unset output
 
 
 
set tit  "Phase advance [2pi] versus s "  font "roman,25"
set xlab "s (m)"         font "roman,25"
set ylab "mu_x/2pi, mu_y/2pi (m)" font "roman,25"
plot \
 "betaFromMatrix.out"       us ($1/100.):($10) w l lw 2  tit "mu_x" , \
 "betaFromMatrix.out"       us ($1/100.):($11) w l lw 2  tit "mu_y" 
 
 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_muxy.eps"
 replot
 set terminal X11
 unset output
 
exit
 
