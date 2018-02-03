
set title "Plotted from file zgoubi.plt  \n From polar frame to lab. (X,Y)"

set key maxcol 1
set key t r

#set logscale y 

set xtics mirror
set ytics mirror

set xlabel 'Y * cos(X)  [m]'
set ylabel 'Y * sin(X)  [m]'

cm2m = 0.01
MeV2eV = 1e6
am = 938.27203
c = 2.99792458e8

set xrange []
set x2range []

plot  \
   'zgoubi.plt' u ($10 *cm2m *cos($22)):($10 *cm2m *sin($22)) w l lc rgb 'red' tit 'B vs. x_lab, y_lab'

     set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_zgoubi.plt_XYLab.eps"  
       replot  
       set terminal X11  
       unset output  

 
pause 1

plot  \
   'zgoubi.plt' u ($22):($10) w l lc rgb 'red' tit 'R vs. angle'

pause 88
exit

