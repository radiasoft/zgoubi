
 set xlabel "x (cm)" font "roman,18"
 set ylabel "y (cm)" font "roman,18"
# set title "Dynamic aperture" font "roman,20"

 set xtics font "roman,14"
 set ytics font "roman,14"


plot  \
  "searchDA.out_save"  using 1:($4==1 ? $2 : 1/0) w l lw 2 lt 1 title 'dp/p= -0.2%'    , \
  "searchDA.out_save"  using 1:($4==4 ? $2 : 1/0) w l lw 2 lt 2 title 'dp/p= -0.1%'    , \
  "searchDA.out_save"  using 1:($4==2 ? $2 : 1/0) w l lw 2 lt 3 title 'dp/p= 0', \
  "searchDA.out_save"  using 1:($4==5 ? $2 : 1/0) w l lw 2 lt 4 title 'dp/p= +0.1%'    , \
  "searchDA.out_save"  using 1:($4==3 ? $2 : 1/0) w l lw 2 lt 5 title 'dp/p= +0.2%'

pause 8


 set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12
 set output "gnuplot_DAs.eps"
 replot
 set terminal X11
 unset output

