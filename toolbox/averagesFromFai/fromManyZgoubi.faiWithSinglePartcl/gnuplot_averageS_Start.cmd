
 set xlabel "G.\gamma"       font "roman,20"
 set ylabel "<S_y>  "               font "roman,20"
 set title "Zgoubi tracking (Nov. 2011) " font "roman,18"
 
       set bmargin 4.5
       set lmargin 10
       set rmargin 4
       set tmargin 3

 set xtics font "roman,14"
 set ytics font "roman,14"

#set xrange [43.5:49.1]   # AGS
#set xrange [:394]   # RHIC

avSy(x) = a*x + avSy0
fit_limit=1e-12

avSy(x) = a*x**2 + b*x + c
fit_limit=1e-12
#######  231+Qy
#fit avSy(x) \
#     'averageS.Out' u ($8):($9==2000 && $2>0 && $8>256 && $8<257.8 ? $6 : 1/0) via a, b, c
#plot \
#     'averageS.Out' u ($8):($9==2000 && $2>0 && $8<258 ? $6 : 1/0) w l lw 1. tit '4000 prtcls' , \
#      avSy(x) w l lw 1.5 lt 3 tit "fit, Average Sy, Final"
######  411-Qy
#fit avSy(x) \
#     'averageS.Out' u ($8):($9==2000 && $2>0 && $8>375 && $8<377 ? $6 : 1/0) via a, b, c
#plot \
#     'averageS.Out' u ($8):($9==2000 && $2>0 &&  $8<377 ? $6 : 1/0) w l lw 1. tit '4000 prtcls' , \
#     avSy(x) w l lw 1.5 lt 3 tit "fit, Average Sy, Final"
#######  393+Qy
fit avSy(x) \
     'averageS.Out' u ($8):($9==2000 && $2>0 && $2>0 && $8<415.2 ? $6 : 1/0) via a, b, c
plot \
     'averageS.Out' u ($8):($9==2000 && $2>0 && $2>0 && $8<416. ? $6 : 1/0) w l lw 1. tit '4000 prtcls' , \
     avSy(x) w l lw 1.5 lt 3 tit "fit, Average Sy, Final"

 set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12
 set output "gnuplot_averageSy_Start.eps"
 replot
 set terminal X11
 unset output



print ' '
print 'minimum : ', avSy(-b/2/a), '   location, G.gamma= ',-b/2/a


pause 10


exit


