
 set xlabel "x (cm)" font "roman,18"
 set ylabel "y (cm)" font "roman,18"
 set title "Dynamic aperture" font "roman,20"

 set xtics font "roman,14"
 set ytics font "roman,14"


plot  \
  "searchDA.out"  using 1:($4==1 ? $2 : 1/0) w l lw 2 lt 1 title '10 MeV', \
  "searchDA.out"  using 1:($4==2 ? $2 : 1/0) w l lw 2 lt 2 title '12 MeV', \
  "searchDA.out"  using 1:($4==3 ? $2 : 1/0) w l lw 2 lt 3 title ' ', \
  "searchDA.out"  using 1:($4==4 ? $2 : 1/0) w l lw 4 lt 4 title ' ', \
  "searchDA.out"  using 1:($4==5 ? $2 : 1/0) w l lw 2 lt 5 title ' ', \
  "searchDA.out"  using 1:($4==6 ? $2 : 1/0) w l lw 2 lt 6 title ' ', \
  "searchDA.out"  using 1:($4==7 ? $2 : 1/0) w l lw 2 lt 7 title ' '

pause 20

 set terminal postscript eps blacktext color
 set output "gnuplot_DAs.eps"
#set output '| cat >> gnuplot_DAs.eps'
 replot
 set terminal X11
 set output
