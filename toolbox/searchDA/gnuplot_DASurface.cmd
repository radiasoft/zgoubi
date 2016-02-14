
set title "1000-cell admittance vs. momentum \n using QF, BD OPERA maps"

 set xlabel "p/p_0 (p_0=67 MeV)" font "roman,18"
 set ylabel "Admittance (mm^2)" font "roman,18"

 set xtics font "roman,14"
 set ytics font "roman,14"

cmsq2mmsq=1.e2

plot   \
   'surface.out'  using ($1):($2*cmsq2mmsq) w lp pt 5 ps .4 lw 2 lt 1 lc 3  notit    

pause 8


 set terminal postscript eps blacktext color enh size 8.3cm,5cm "Times-Roman" 12
 set output "gnuplot_DASurface.eps"
 replot
 set terminal X11
 unset output

