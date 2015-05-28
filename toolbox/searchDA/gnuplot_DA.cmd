
set tit "C{/Symbol b} cell, 1000-cell DA  \n (QF-Displaced_{\261}QD-CFM.table field map) \n E1-5 : 76, 146, 216, 225, 286 MeV" font "roman,14 

 set xlabel "y (mm)" font "roman,14"
 set ylabel "z (mm)" font "roman,14"

 set xtics font "roman,12"
 set ytics font "roman,12"

# set bmargin -1

cm2mm=10.

#set xrange [-12: 12]

plot for [i=1:5] \
  'searchDA.out_save'  using (($1+$5)*cm2mm):($4==i ? $2*cm2mm : 1/0) w lp pt i ps .4 lw 2 lt 1 lc i title 'E'.i 


 set terminal postscript eps blacktext color enh size 8.3cm,5cm "Times-Roman" 12
 set output "gnuplot_DAs.eps"
 replot
 set terminal X11
 unset output

pause 8


