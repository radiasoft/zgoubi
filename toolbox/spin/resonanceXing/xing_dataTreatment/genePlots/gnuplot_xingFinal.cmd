 set xlabel "Turn"       font "roman,20"
 set ylabel "Sy "               font "roman,20"
 set title "RESONANCE CROSSING, LAST ~200 TURNS" font "roman,18"

 set xtics font "roman,14"
 set ytics font "roman,14"

  plot 'fromBFai2Fai.out' us ($38):($38>maxPnt && $38 < maxPnt+30 ?   $22 : 1/0)  w l lt 3 lw .5  title 'Zgoubi' 

#pause .1 

 set terminal postscript eps blacktext color
 set output "gnuplot_xingFinal.eps"
 replot
 set terminal X11
 unset output

exit





plot  'tunesFromMatrix.out_HER' us ($10):($3>0 ? $3 : 1/0) w l lw 2 lt 1 tit 'Qx-45' , \
      'tunesFromMatrix.out_HER' us ($10):($3>0 ? $3 : 1/0) w p ps 1 pt 5 tit 'Qx-45' , \
      'tunesFromMatrix.out_HER' us ($10):($4>0 ? $4 : 1/0) w l lw 2 lt 2 tit 'Qy-20' , \
      'tunesFromMatrix.out_HER' us ($10):($4>0 ? $4 : 1/0) w p ps 1 pt 6 tit 'Qy-20'

#plot [] [-0.2:0.2] 'tunesFromMatrix.out' us ($10):($3>0 ? ($3-0.538990) : 1/0) w l lw 2 lt 1 tit 'Qx-Qx_0' , \
#                   'tunesFromMatrix.out' us ($10):($4>0 ? ($4-0.561648) : 1/0) w l lw 2 lt 2 tit 'Qy-Qy_0'

pause 5

 set terminal postscript eps blacktext color
 set output "gnuplot_TunesFromMatrix.eps"
 replot
 set terminal X11
 set output



G = 1.7928474
qz = 8.7746
Jn2=0.336e-5

# x = gG - gGRso = gG - (4M - qz)      gGRso=4M-qz ~ 39

rho2(x) = 1. -( ( (x-(48-qz) )*(x-(48-qz)) / ( (x-(48-qz))*(x-(48-qz))  + Jn2 ) ))
# plot [39:40] S(x)

qz = 8.7746
Jn2=0.336e-5

fit_limit=1e-12

fit rho2(x) 'gnu_fitStatic.data' using (G*$10):(1.-($5+$6)*($5+$6)/4) via 'gnu_fitStart'

set y2tics
set y2range [0.05:2.5]

 plot  [-0.025:0.025] 'gnu_fitStatic.data' us (G*$10 -(48 -qz)):((1.-($5+$6)*($5+$6)/4)) smooth cspline w l tit 'Smooth', \
(rho2(x+(48-qz))) w l lt 1 lw 2 tit 'rho', \
'gnu_fitStatic.data' us (G*$10 -(48 -qz)):((1.-($5+$6)*($5+$6)/4))  w p ps 2 tit 'Zgoubi', \


'gnu_fitStatic.data' us (G*$10 -(48 -qz)):((1.-($5+$6)*($5+$6)/4))  w p ps 2 tit 'Zgoubi', \
'gnu_fitStatic.data' us (G*$10 -(48 -qz)):(1e5*(G*$10-(48 -qz))**2/(1/(1-($5+$6)*($5+$6)/4)-1)) axes x1y2 w l lt 3 lw .5 tit '|J_n| num.'

pause 2

 set terminal postscript eps blacktext color
 set output "gnuplot_fitStatic_rho.eps"
 replot
 set terminal X11
 set output

plot  1e5*(x-(48-qz))**2/(1./rho2(x) -1.)  w l lt 4 lw .5 tit '|Jn^2| / fit'


