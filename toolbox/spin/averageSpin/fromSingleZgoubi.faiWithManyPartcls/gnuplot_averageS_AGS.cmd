
 set xlabel "G.gamma"       font "roman,14"
 set x2label "G.gamma"       font "roman,14"
 set ylabel "<S_y>  "               font "roman,14"
 set y2label "<S_y> in Gg=10.5 region "               font "roman,14"
 set title "Zgoubi tracking (Feb. 2012) " font "roman,14"
 
       set bmargin 4.5
       set lmargin 10
       set rmargin 15
       set tmargin 3

 set xtics font "roman,14" nomirror
 set x2tics font "roman,10" nomirror
 set ytics font "roman,14" nomirror
 set y2tics font "roman,10" nomirror

#set xrange [4.5:10.]   # AGS
#set xrange [:394]   # RHIC

plot \
     'averageS.out' u ($1):($2) axes x2y2  w p pt 2 ps .2 notit  
#     'averageS.out' u ($1):($3>0 && $4<25000 ? $2 : 1/0) axes x1y1 w l lw .8 notit 
#     'averageS.out' u ($1):($1>10.4 && $1<10.6 ? $2 : 1/0) axes x2y2  w p pt 2 ps .2 notit  

 set terminal postscript eps blacktext color enh size 8.3cm,4.4cm "Times-Roman" 14
 set output "gnuplot_averageSy.eps"
 replot
 set terminal X11
 unset output

pause 10


exit


