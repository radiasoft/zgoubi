 set xtics font "roman,12"
 set ytics font "roman,12"
 
set tit  "beta_{r,y} versus p/p_0" font "roman,20"
set xlab "p/p_0"             font "roman,16"
set ylab "beta_r, beta_y (m)"  font "roman,16"
plot \
 "betaFromMatrix.out" us ($14):($3) w l lw 2  tit "beta_r" , \
 "betaFromMatrix.out" us ($14):($7 * (1)) w l lw 2  tit "beta_y" 
 
 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_betxy.eps"
 replot
 set terminal X11
 unset output
 
 
set tit  "D_r, D_y versus p/p_0" font "roman,20"
set xlab "p/p_0"             font "roman,16"
set ylab "D_r, D_y (m)"      font "roman,16"
plot \
 "betaFromMatrix.out" us ($14):($4) w l lw 2  tit "D_r" , \
 "betaFromMatrix.out" us ($14):($8) w l lw 2  tit "D_y" 
 
 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_Dxy.eps"
 replot
 set terminal X11
 unset output
 
 
 
set tit  "Fractional tunes versus p/p_0 "  font "roman,20"
set xlab "p/p_0"         font "roman,16"
set ylab "Q_r, Q_y" font "roman,16"
plot \
 "betaFromMatrix.out"       us ($14):($10) w l lw 2  tit "Q_r" , \
 "betaFromMatrix.out"       us ($14):($11) w l lw 2  tit "Q_y" 
 
 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_muxy.eps"
 replot
 set terminal X11
 unset output

pause 88 
exit
 
