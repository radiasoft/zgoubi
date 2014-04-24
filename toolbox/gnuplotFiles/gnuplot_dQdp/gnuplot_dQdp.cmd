 
Qx0 = 0.196091
Qy0 = 0.192599
L0 = 8.12802050E+04
T0 = 2.71121590

dEE = 1.e-3 + 0.000000001   #just to get the horizontal scale span correct
set xlabel "dE/E / 10^{-3}" font "roman,18"
set title "Chromaticity" font "roman,20"

 set xtics font "roman,12"
 set ytics font "roman,12"

 set ylabel "dQ" font "roman,18"
plot "tunesFromMatrix.out" u (($13-1.)/dEE):(($3-Qx0))     w l lt 1 lw 3 tit 'dQx', \
     "tunesFromMatrix.out" u (($13-1.)/dEE):(($4-Qy0)) w l lt 3 lw 1 tit 'dQy'

 set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12
 set output "gnuplot_chroma.eps"
 replot
 set terminal X11
 unset output

dQ(x) = dQ0 + dQ1 * x + dQ2 * x*x + dQ3 * x*x*x

fit_limit=1e-12

fit dQ(x) \
     "tunesFromMatrix.out" u (($13-1.)):(($3-Qx0))     via 'gnu_fitStart'

pause 3

fit dQ(x) \
     "tunesFromMatrix.out" u (($13-1.)):(($4-Qy0)) via 'gnu_fitStart'

#-------------- 
set title "Momentum compaction" font "roman,20"

 set ylabel "alpha * dp/p" font "roman,18"
plot "tunesFromMatrix.out" u (($13-1.)/dEE):((($11-L0)/L0))     w l lt 1 lw 1 tit 'alpha'

pause 3 

 set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12
 set output "gnuplot_alpha.eps"
 replot
 set terminal X11
 unset output

alp(x) = alp0 + alp1 * x + alp2 * x*x + alp3 * x*x*x

fit_limit=1e-12

fit alp(x) \
     "tunesFromMatrix.out" u (($13-1.)):(($11-L0)/L0)     via 'gnu_fitStart_alp'

#-------------- 
set title "Phase slip factor" font "roman,20"

 set ylabel "eta * dp/p" font "roman,18"
plot "tunesFromMatrix.out" u (($13-1.)/dEE):((($10-T0)/T0))     w l lt 1 lw 1 tit 'eta'

pause 3 

 set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12
 set output "gnuplot_eta.eps"
 replot
 set terminal X11
 unset output

eta(x) = eta0 + eta1 * x + eta2 * x*x + eta3 * x*x*x


fit_limit=1e-12

fit eta(x) \
     "tunesFromMatrix.out"  u (($13-1.)/dEE):((($10-T0)/T0))     via 'gnu_fitStart_eta'


