 
Qx0 = 0.196091
Qy0 = 0.192599

dEE = 1.e-3 + 0.000000001   #just to get the horizontal scale span correct
set xlabel " / sigma" font "roman,18"
set title "Amplitude detuning" font "roman,20"

 set xtics font "roman,12"
 set ytics font "roman,12"

 set ylabel "dQ" font "roman,18"
plot "tunesFromFai.out" u (($13-1.)/dEE):(($3-Qx0))     w l lt 1 lw 3 tit 'dQx', \
     "tunesFromFai.out" u (($13-1.)/dEE):(($4-Qy0)) w l lt 3 lw 1 tit 'dQy'



 set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12
 set output "gnuplot_chroma.eps"
 replot
 set terminal X11
 unset output



dQ(x) = dQ0 + dQ1 * x + dQ2 * x*x + dQ3 * x*x*x

fit_limit=1e-12

fit dQ(x) \
     "tunesFromFai.out" u (($13-1.)):(($3-Qx0))     via 'gnu_fitStart'

pause 3

fit dQ(x) \
     "tunesFromFai.out" u (($13-1.)):(($4-Qy0)) via 'gnu_fitStart'


exit 
