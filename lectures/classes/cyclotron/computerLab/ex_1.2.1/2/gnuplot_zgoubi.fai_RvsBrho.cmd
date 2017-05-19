


set key maxcol 1
set key t l

#set logscale y 

set tit 'R vs. particle rigidity, from zgoubi.fai'

set xtics nomirror
set x2tics nomirror
set ytics 

set xlabel 'B{/Symbol r}   /T.m'
set x2label 'kin. E   /MeV'
set ylabel 'R   /cm'

### electron : 
#am = 0.511
#### proton :
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

set xrange [Brmi:Brma]
set x2range [Ekmi:Ekma]

plot  \
   'zgoubi.fai' u (sqrt($25**2-$29**2)*MeV2eV/c):10 w p pt 5 ps .8 lc 1  tit "R(B{/Symbol r})",\
    rho(x) w l lc 1 tit 'R(B{/Symbol r}), theor.' ,\
   'zgoubi.fai' u ($25-$29):10 axes x2y1 w p pt 5 ps 0.4 lc 3 tit "R(E_kin)"

     set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_zgoubi.fai_TOFvsBrho.eps"  
       replot  
       set terminal X11  
       unset output  

 
pause 1
exit

