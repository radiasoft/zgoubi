


set key maxcol 1
set key t l

#set logscale y 

set tit 'BZ vs. particle s, from zgoubi.plt'

set xtics 
set ytics 

set xlabel 's   /cm'
set ylabel 'BZ   /kG'

#! electron : 
#am = 0.511
#! proton :
am = 938.27203
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

#set xrange [Brmi:Brma]
#set x2range [Ekmi:Ekma]

plot  \
   'zgoubi.plt' u ($14):($25) w p pt 5 ps .4 lc 1  notit 

     set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_zgoubi.plt_BZ.vs.s.eps"  
       replot  
       set terminal X11  
       unset output  

 
pause 2
exit

