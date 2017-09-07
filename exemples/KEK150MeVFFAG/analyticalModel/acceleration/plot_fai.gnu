
set tit 'Kinetic energy versus RF phase'
set xlab 'phase (radian)'
set ylab 'E_K [MeV]'

# header lines in zgoubi.fai tell that phase (RET, for the French 'retard') 
# and kinetic-energy are respectivewly  columns 24 and 34  

plot  'zgoubi.fai' us 34:24 w p ps .5 tit "longitudinal space"

 set terminal postscript eps blacktext color
 set output "plot_fai_phase-Ek.eps"
 replot
 set terminal X11
 set output

pause 2

set tit 'Horizontal phase-space'
set xlab 'x (cm)'
set ylab "x' (rad)"

# header lines in zgoubi.fai tell that x (Y) and x' (T)  are respectivewly  columns 2 and 3

plot  'zgoubi.fai' us 12:13 w p ps .8 tit "damped xx' space"

 set terminal postscript eps blacktext color
 set output "plot_fai_YT.eps"
 replot
 set terminal X11
 set output

pause 2

exit
