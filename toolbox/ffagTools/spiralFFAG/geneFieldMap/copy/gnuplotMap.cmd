
set xlabel 'X'
set ylabel 'Y'
set zlabel 'B (arb. units)'
set title "Field map of spiral dipole"

# ticslevel : (pos-zmin)/(zmin-zmax) is for adjusting the z offset of the plot
#set ticslevel (0. -0.)/(0. - 2.)
 
#set hidden3d
set view 39, 257
# splot "fort.88"  with lines 
 splot "fort.88"  with dots 

pause 1

 set terminal postscript eps blacktext color
 set output "gnuplot.eps"
 replot

set terminal X11
set output

#unset hidden3d
unset view 


