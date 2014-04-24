
#set title 'TUNE DIAGRAM' font "helvetica,16"
set xlab 'Qx' font "helvetica,22"
set ylab 'Qy' font "helvetica,22"
plot  \
      'tunesFromMatrix' us ($3):($4) w l lw 2 lt 1 tit 'Hard edge' , \
      'gnuplotTuneDiagWithLines.data' us ($1):($2) w l lw .5 lt 3 notit 

pause 2

set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set output "gnuplot_QxQy.eps"
 replot
 set terminal X11
 unset output

exit
