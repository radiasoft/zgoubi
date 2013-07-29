reset
 
set xlabel "Turn #" font "roman,18"
set ylabel "Emittances, x, y" font "roman,18"
set y2label "\epsilon_y/\epsilon_x" font "roman,18" 
set title "EMITTANCES" font "roman,20"
 
set xtics nomirror font "roman,11"
set ytics nomirror font "roman,11"
set y2tics nomirror font "roman,11"
 
set key font "roman,10" center right

set logscale y 10 
set logscale y2 10 

  FIT_LIMIT = 1e-6

ex(x) = exf + (ex0-exf) * exp(-ax*x)
ax               = 1e-3
ex0              = 1e-8
exf              = 1e-8

ey(x) = eyf + (ey0-eyf) * exp(-ay*x)
ay               = 1.e-3
ey0              = 1e-8
eyf              = 1e-12

eyx(x) = eyxf
eyxf              = 0.1

set xrange [15000:20000]
fit eyx(x)  "./lipsFitFromFai_iterate.Out"  u 11:($4/$3) via eyxf
set xrange [110:20000]
fit ex(x)   "./lipsFitFromFai_iterate.Out"  u 11:3 via ax, ex0, exf
set xrange [110:20000]
fit ey(x)   "./lipsFitFromFai_iterate.Out"  u 11:4 via ay, ey0, eyf

plot \
    "./lipsFitFromFai_iterate.Out"    u 11:3 w p ps .8 tit 'E_x' ,\
    ex(x) w l lw 1.9 tit 'E_x, fit' ,\
    "./lipsFitFromFai_iterate.Out"    u 11:4 w p ps .8 tit 'E_y' ,\
    ey(x) w l lw 1.9 tit 'E_y, fit' ,\
    "./lipsFitFromFai_iterate.Out"    u 11:($4/$3) axes x1y2 w p ps .8  tit 'E_y/E_x' ,\
      eyx(x) axes x1y2 w l lw 1.9 tit 'E_y/E_x, fit' 

set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_emittances.eps" 
 replot 
 set terminal X11 
 unset output 

print " tau_x=  ",1/ax , "   ex_f=  ",exf,  "  tau_y=  ",1/ay, "   ey_f=  ", eyf, "   ey_f/ex_f=  ",eyxf, "  ex0=  ", ex0

pause 888
exit
##########################################3
##########################################3
##########################################3

