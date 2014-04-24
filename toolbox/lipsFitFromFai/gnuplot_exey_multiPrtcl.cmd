# ex/ex0 and ey/ex
set xlabel "C/Delta" font "roman,18"
set ylabel "Emittance ratio" font "roman,18"
set title "EMITTANCE RATIOS" font "roman,20"
 
set xtics font "roman,12"
set ytics font "roman,12"
 
set key font  "roman,12"  center left
set logscale x  10

Da = 0.01
Db = 400.
Dc = 0.01

#ex0 = 6.79e-09
ex0 = 6.66e-09

fac = 1.

ee0(x) = (2.* (x/fac/Da)*(x/fac/Da)+ 1.) / (4. *(x/fac/Da)*(x/fac/Da) + 1.) 
eyex(x) = (x/fac/Db)*(x/fac/Db) / ((x/fac/Db)*(x/fac/Db) + 0.5) + Dc

fit_limit=1e-12
fit ee0(x) \
   "./emittances.out" u (fac*$1):($5)  via  Da
fit eyex(x) \
   "./emittances.out" u (fac*$1):($11)  via  Db, Dc
      
#set xrange[0:.52]
#set yrange [-.02:1.02]

plot \
   "./emittances.out" u (fac*$1):($5/ex0) w p ps 1. tit 'ex/ex0',\
   ee0(x) w l lw 1. tit 'e_x/e_x0, fit'  ,\
   "./emittances.out" u (fac*$1):($11) w p ps 1. tit 'ey/ex',\
   eyex(x) w l lw 1. tit 'e_y/e_x, fit' 

print Da, Db

 set term post eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12
 set output "gnuplot_exey_multiPrtcl.eps"
 replot
 set terminal X11
 unset output

pause 8888888
