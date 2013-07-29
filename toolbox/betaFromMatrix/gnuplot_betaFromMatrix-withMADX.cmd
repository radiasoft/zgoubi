 set xtics font "roman,18"
 set ytics font "roman,18"
 
set tit  "beta_x" font "roman,25"
set xlab "s (m)"             font "roman,25"
set ylab "beta_x (m)"  font "roman,25"

set xrange[0:807.1/4]
set yrange[9:24]
plot \
 "betaFromMatrix.out" us ($1/100.):($3) w l lt 1 lc 3 lw 2  tit "Zgoubi" , \
 "twiss4zgoubi" us ($12):($14) w l lt 1 lc 2 lw 2  tit "MADX" 
 
 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_betx_ZMAD.eps"
 replot
 set terminal X11
 unset output
 
 
set tit  "beta_y" font "roman,25"
set xlab "s (m)"    font "roman,25"
set ylab "beta_y (m)"  font "roman,25"

set xrange[0:807.1/4]
set yrange[9:24]
plot \
 "betaFromMatrix.out" us ($1/100.):($7) w l lt 1 lc 3 lw 2  tit "Zgoubi" , \
 "twiss4zgoubi" us ($12):($16) w l lt 1 lc 2 lw 2  tit "MADX" 
 
 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_bety_ZMAD.eps"
 replot
 set terminal X11
 unset output
 
 
pause 1

reset

set tit  "X orbit" font "roman,25"
set xlab "s (m)"             font "roman,25"
set ylab "x (m)"      font "roman,25"

set xrange[0.01:807.]
plot \
 "betaFromMatrix.out" us ($1/100.):($16/100.) w l lt 1 lc 3 lw 2  tit "Zgoubi" , \
 "twiss4zgoubi" us ($12):($17) w l lt 1 lc 2 lw 2  tit "MADX" 
 
 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_x_ZMAD.eps"
 replot
 set terminal X11
 unset output
 
 
set tit  "Y orbit" font "roman,25"
set xlab "s (m)"             font "roman,25"
set ylab "D_y (m)"      font "roman,25"

set xrange[0:807.1]
plot \
 "betaFromMatrix.out" us ($1/100.):($18/100.) w l lt 1 lc 3 lw 2  tit "Zgoubi" , \
 "twiss4zgoubi" us ($12):($18) w l lt 1 lc 2 lw 2  tit "MADX" 
 
 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_y_ZMAD.eps"
 replot
 set terminal X11
 unset output
 
pause 1

reset

set tit  "D_x" font "roman,25"
set xlab "s (m)"             font "roman,25"
set ylab "D_x (m)"      font "roman,25"

set xrange[0:807.1/4]
plot \
 "betaFromMatrix.out" us ($1/100.):($4) w l lt 1 lc 3 lw 2  tit "Zgoubi" , \
 "twiss4zgoubi" us ($12):($19) w l lt 1 lc 2 lw 2  tit "MADX" 
 
 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_Dx_ZMAD.eps"
 replot
 set terminal X11
 unset output
 
 
set tit  "D_y" font "roman,25"
set xlab "s (m)"             font "roman,25"
set ylab "D_y (m)"      font "roman,25"

set xrange[0:807.1/4]
plot \
 "betaFromMatrix.out" us ($1/100.):($8) w l lt 1 lc 3 lw 2  tit "Zgoubi" , \
 "twiss4zgoubi" us ($12):($21) w l lt 1 lc 2 lw 2  tit "MADX" 
 
 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_Dy_ZMAD.eps"
 replot
 set terminal X11
 unset output
 
pause 1

reset
 
set tit  "Phase advance [2pi]"  font "roman,25"
set xlab "s (m)"         font "roman,25"
set ylab "mu_x/2pi" font "roman,25"

set xrange[0:807.1/4]
plot \
 "betaFromMatrix.out"       us ($1/100.):($12) w l lt 1 lc 3 lw 2  tit "Zgoubi" , \
 "twiss4zgoubi" us ($12):($23 - int($23)) w l lt 1 lc 2 lw 2  tit "MADX" 

 
 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_mux_ZMAD.eps"
 replot
 set terminal X11
 unset output

set tit  "Phase advance [2pi]"  font "roman,25"
set xlab "s (m)"         font "roman,25"
set ylab "mu_y/2pi" font "roman,25"

set xrange[0:807.1/4]
plot \
 "betaFromMatrix.out"       us ($1/100.):($13) w l lt 1 lc 3 lw 2  tit "Zgoubi" , \
 "twiss4zgoubi" us ($12):($24 - int($24)) w l lt 1 lc 2 lw 2  tit "MADX" 
 
 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_muy_ZMAD.eps"
 replot
 set terminal X11
 unset output

 
pause 1
 
exit
 
