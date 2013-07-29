# H and V CLOSED ORBITS 
set xlabel "s (m)" font "roman,18"
set xtics font "roman,12"
set ytics font "roman,12"
 
set grid
 
set ylabel "x, y (m)" font "roman,18"
set title "H & V closed orbits along ring" font "roman,20"
 
plot  \
"getDiffFromFai.out"  u ($1 /100.):($2 /100.) w l lw 2  tit "x_{co}", \
"getDiffFromFai.out"  u ($1 /100.):($3 /100.) w l lw 2  tit "y_{co}"
 
pause 3
 
set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_CO.eps" 
 replot 
 set terminal X11 
 unset output 
 
 
set ylabel "D_x, D_y (m)" font "roman,18"
set title "H & V dispersion along ring" font "roman,20"
 
plot  \
"getDiffFromFai.out"  u ($1 /100.):($4 /100.) w l lw 2  tit "D_{x}",  \
 "getDiffFromFai.out"  u ($1 /100.):($5 /100.) w l lw 2  tit "D_{y}"
 
pause 3
 
set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_DxDy.eps" 
 replot 
 set terminal X11 
 unset output 
 
exit 
