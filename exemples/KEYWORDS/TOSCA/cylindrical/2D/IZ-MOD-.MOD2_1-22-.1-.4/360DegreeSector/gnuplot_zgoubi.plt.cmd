
set title "Plotted from file zgoubi.plt  \n From zgoubi's polar frame to lab frame \n u ($10 *cm2m *cos($22)):($10 *cm2m *sin($22)) " font "sans, 14"

set key maxcol 1
set key c c 

#set logscale y 

set xtics mirror font  "sans, 14"
set ytics mirror font  "sans, 14"

set size ratio 1

set xlabel 'Y * cos(X)  [m]' font  "sans, 14"
set ylabel 'Y * sin(X)  [m]' font  "sans, 14"

cm2m = 0.01
MeV2eV = 1e6
am = 938.27203
c = 2.99792458e8

set xrange []
set x2range []

plot  \
   'zgoubi.plt' u ($10 *cm2m *cos($22)):($10 *cm2m *sin($22)) w l lc rgb 'red' tit 'x\_lab, y\_lab'

#     set terminal postscript eps blacktext color  enh size 8cm,8cm "Times-Sans" 12  
     set terminal postscript eps blacktext color  enh  "Times-Sans" 12  
       set output "gnuplot_zgoubi.plt_XYLab.eps"  
       replot  
       set terminal X11  
       unset output  

 
pause 3
exit

