
set grid
       set bmargin 4
       set lmargin 11
       set rmargin 11
       set tmargin 3

 set key center top font "roman,20"

 set xtics font "roman,18"
 set ytics font "roman,18"
 
set tit  "beta_r,z versus p/p_{inj}" font "roman,20"
set xlab "p/p_{inj}"             font "roman,20"
set ylab "beta_r, beta_z (m)"  font "roman,20"
plot  \
 "betaFromMatrix.out" us ($14 >0 ? $14 : 1/0):($3) w l lw 2.5  tit "\beta_r" , \
 "betaFromMatrix.out" us ($14 >0 ? $14 : 1/0):($7) w l lw 2.5  tit "\beta_z" 
 
 set term post eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_betrz.eps"
 replot
 set terminal X11
 unset output
 
 
set tit  "D_r, D_z versus p/p_{inj}" font "roman,20"
set xlab "p/p_{inj}"             font "roman,20"
set ylab "D_r, D_z (m)"      font "roman,20"
plot  \
 "betaFromMatrix.out" us ($14 >0 ? $14 : 1/0):($4) w l lw 2.5  tit "D_r" , \
 "betaFromMatrix.out" us ($14 >0 ? $14 : 1/0):($8) w l lw 2.5  tit "D_z" 
 
 set term post eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_Drz.eps"
 replot
 set terminal X11
 unset output
 
 
 
set tit  "Tunes versus p/p_{inj}"  font "roman,20"
set xlab "p/p_{inj}"         font "roman,20"
set ylab "Qr,  Qz" font "roman,20"
plot \
 "betaFromMatrix.out"       us ($14 >0 ? $14 : 1/0):($12) w l lw 2.5  tit "Qr" , \
 "betaFromMatrix.out"       us ($14 >0 ? $14 : 1/0):($13 < .4 ? $13 : $13-1) w l lw 2.5  tit "Qz" 
 
 set term post eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_murz.eps"
 replot
 set terminal X11
 unset output
 
exit
 
