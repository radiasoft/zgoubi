
set grid

set xlabel "Kin. E (MeV)" font"roman, 22"
set ylabel "Qx, Qy"  font"roman, 22"
 
set key top center

am = 0.511
am2 = am*am

plot   [][] \
     "tunesFromFai.out"  u (sqrt($17*$17 + am2)-am):($5) w p pt 4 ps .4 lw 2 lc 1  tit "Qx" ,\
     "tunesFromFai.out"  u (sqrt($17*$17 + am2)-am):($6) w p pt 4 ps .4 lw 2 lc 3  tit "Qy" 

 set terminal postscript eps blacktext color
 set output "gnuplot_tunes.eps"
 replot
 set term X11

 print " " 
 print "Tune scan graphic saved in ./gnuplot_tunes.eps"
 print " " 

pause 2. 

unset grid
set xlabel "Q_x" font"roman, 22"
set ylabel "Q_y" font"roman, 22"

set key maxcol 1
set key t l

plot  [][] \
     "tunesFromFai.out"  u ($12>= 1 && $12 <= 5 ? $5 : 1/0):($6) w p pt  4 ps 1.5 lw 3 lc 1  tit  " 40 MeV"  ,\
     "tunesFromFai.out"  u ($12>= 6 && $12 <=10 ? $5 : 1/0):($6) w p pt  6 ps 1.5 lw 3 lc 2  tit  " 42 MeV"  ,\
     "tunesFromFai.out"  u ($12>=11 && $12 <=15 ? $5 : 1/0):($6) w p pt  8 ps 1.5 lw 3 lc 3  tit  " 78 MeV"  ,\
     "tunesFromFai.out"  u ($12>=16 && $12 <=20 ? $5 : 1/0):($6) w p pt 10 ps 1.5 lw 3 lc 4  tit  "114 MeV"  ,\
     "tunesFromFai.out"  u ($12>=21 && $12 <=25 ? $5 : 1/0):($6) w p pt 12 ps 1.5 lw 3 lc 5  tit  "150 MeV"  ,\
     "tunesFromFai.out"  u ($12>=26 && $12 <=30 ? $5 : 1/0):($6) w p pt 14 ps 1.5 lw 3 lc 6  tit  "160 MeV"  ,\
     "tunesDiagLines_x0-.5_y0-.5.data" u 1:2 w l lc 3 notit

 set terminal postscript eps blacktext color
 set output "gnuplot_tuneDiag.eps"
 replot
 set term X11

 print " " 
 print "Tune diagram graphic saved in ./gnuplot_tuneDiag.eps"
 print " " 

pause 4.

exit




