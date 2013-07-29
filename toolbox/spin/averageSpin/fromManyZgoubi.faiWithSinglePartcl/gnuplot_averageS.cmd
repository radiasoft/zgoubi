
 set xlabel "G.gamma"       font "roman,20"
 set ylabel "<S_y>  "               font "roman,20"
 set title "Zgoubi tracking (May 2012) " font "roman,18"
 
       set bmargin 4.5
       set lmargin 10
       set rmargin 4
       set tmargin 3

 set xtics font "roman,14"
 set ytics font "roman,14"

#set xrange [43.5:49.1]   # AGS
#set xrange [:394]   # RHIC

plot \
     'averageS.Out' u ($8):($9==1290 && $2>0 ? $6 : 1/0) w l lw .5 tit 'no cutoff 2k prtcls'  
#     'averageS.Out_full' u ($8):($9==4032 && $2>0 ? $6 : 1/0) w l lw .5 tit 'no cutoff 2k prtcls'  

 set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12
 set output "gnuplot_averageSy.eps"
 replot
 set terminal X11
 unset output

pause 888


exit


