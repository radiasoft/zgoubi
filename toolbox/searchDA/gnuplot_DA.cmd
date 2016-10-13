
set title "1000-cell admittance \n using QF, BD OPERA maps"

 set xlabel "x (mm)" font "roman,18"
 set ylabel "y (mm)" font "roman,18"

 set xtics font "roman,14"
 set ytics font "roman,14"

cm2mm=10.

set xrange [-30:30]

plot   \
   for [i=1:99] 'searchDA.out_save_99trj'  using ($1*cm2mm):($4==i ? $2*cm2mm : 1/0) w lp pt i ps .4 lw 2 lt 1 lc i title 'dp/p='.(i-5).'e-3'    

 set terminal postscript eps blacktext color enh size 8.3cm,5cm "Times-Roman" 12
 set output "gnuplot_DAs.eps"
 replot
 set terminal X11
 unset output

pause 8




