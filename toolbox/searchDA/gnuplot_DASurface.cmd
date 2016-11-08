
set title "1000-cell admittance vs. momentum, at center of 12 cm drift \n Using QF and BD OPERA maps / no defects "

 set xlabel "E (MeV)" font "roman,18"
 set ylabel "Admittance (mm^2)" font "roman,18"

 set xtics font "roman,14"
 set ytics nomirror font "roman,14"
 set y2tics nomirror font "roman,0.1"

set logscale y

set key maxrow 1
set key c b

cmsq2mmsq=1.e2
p0 = 119

set xrange [35:170]
set y2range [:1]

plot   \
   'surface.out'  using ($1*p0):($2*cmsq2mmsq) w lp pt 5 ps .4 lw 2 lt 1 lc 3  tit 'Surf.' ,\
   'designE.data' u ($1):2 axes x1y2 w i lt 1 lc 1 tit 'design energies'

 set terminal postscript eps blacktext color enh size 8.3cm,5cm "Times-Roman" 12
 set output "gnuplot_DASurface.eps"
 replot
 set terminal X11
 unset output

pause 8


