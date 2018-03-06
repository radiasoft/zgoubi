
set key maxcol 1
set key t l

#set logscale y 

set tit 'D=p/p0 vs. turn number, from zgoubi.fai'

set xtics mirror
set ytics nomirror
set y2tics 

set xlabel 'turn'
set ylabel 'D = p/p0-1'
set y2label 'Y (cm)'

! electron : 
am = 0.511
#! proton :
#am = 938.27203

m2cm = 100.
MeV2eV = 1e6
c = 2.99792458e8
B = 0.5   # [T]
rho(x) = m2cm * x/B    # [cm]

Ekmi = 1
Ekma = 5
Emi = Ekmi + am
Ema = Ekma + am
Brmi = sqrt(Emi**2-am**2)*MeV2eV/c
Brma = sqrt(Ema**2-am**2)*MeV2eV/c

set xrange [990:1010]
#set x2range [Ekmi:Ekma]

plot  \
   "turn3000.fai"    u ($38):($9 +1.)     w p pt 4 ps .9 lc rgb "red"  tit "D, track"  ,\
   "zgoubi.fai" u ($38+999):($9 +1.) w p pt 5 ps .6 lc rgb "blue"  tit "D, recvrd" ,\
   "turn3000.fai"    u ($38):($10)     axes x1y2 w p pt 4 ps .9 lc rgb "green"  tit "Y, track"  ,\
   "zgoubi.fai" u ($38+999):($10) axes x1y2 w p pt 5 ps .6 lc rgb "cyan"  tit "recvrd"  

     set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_zgoubi.fai_Turn-D.eps"  
       replot  
       set terminal X11  
       unset output  

 
pause 1
exit

