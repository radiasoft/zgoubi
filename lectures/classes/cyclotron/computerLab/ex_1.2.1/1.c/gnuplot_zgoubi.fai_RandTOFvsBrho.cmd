


set key font "roman, 6"
set key maxrow 4 sample .1 
set key t c

#set logscale y 

set tit 'R and TOF vs. particle rigidity, changing step size'

set xtics nomirror
set x2tics nomirror
set ytics nomirror
set y2tics nomirror

set xlabel 'B{/Symbol r}   /T.m'
set x2label 'kin. E   /MeV'
set ylabel '(R-R_0)/R_0'
set y2label 'TOF   /{/Symbol m}s'

### electron : 
#am = 0.511
#### proton :
am = 938.27203
m2cm = 100.
MeV2eV = 1e6
s2mus = 1e6
c = 2.99792458e8
B = 0.5   # [T]
pi = 4. * atan(1.)
Nsector = 6.

rho(x) = m2cm * x/B    # [cm]
tof(x) = (2. *pi *x/B/Nsector) / ((x*c)/sqrt((x*c)**2 + (am*MeV2eV)**2)) /c *s2mus   # [cm]

Ekmi = 1
Ekma = 5
Emi = Ekmi + am
Ema = Ekma + am
Brmi = sqrt(Emi**2-am**2)*MeV2eV/c
Brma = sqrt(Ema**2-am**2)*MeV2eV/c

set xrange [Brmi:Brma]
set x2range [Ekmi:Ekma]

plot  \
   'zgoubi.fai_stepSize13cm'   u (sqrt($25**2-$29**2)*MeV2eV/c):15 axes x1y2 w lp pt 1  ps .4 lc 1  tit "TOF, step 13 cm" ,\
   'zgoubi.fai_stepSize5cm'    u (sqrt($25**2-$29**2)*MeV2eV/c):15 axes x1y2 w lp pt 2  ps .4 lc 2  tit "step 5 cm" ,\
   'zgoubi.fai_stepSize0.1cm'  u (sqrt($25**2-$29**2)*MeV2eV/c):15 axes x1y2 w lp pt 3  ps .4 lc 3  tit "step 0.1 cm" ,\
   'zgoubi.fai_stepSize0.01cm' u (sqrt($25**2-$29**2)*MeV2eV/c):15 axes x1y2 w lp pt 5  ps .4 lc 4  tit "step 0.01 cm" ,\
   'zgoubi.fai_stepSize13cm'   u (sqrt($25**2-$29**2)*MeV2eV/c):(($10 - $3)/$3)     w lp pt 10 ps .4 lc 1  tit "{/Symbol D}R/R " ,\
   'zgoubi.fai_stepSize5cm'    u (sqrt($25**2-$29**2)*MeV2eV/c):(($10 - $3)/$3)     w lp pt 20 ps .4 lc 2  tit "  " ,\
   'zgoubi.fai_stepSize0.1cm'  u (sqrt($25**2-$29**2)*MeV2eV/c):(($10 - $3)/$3)     w lp pt 30 ps .4 lc 3  tit "  " ,\
   'zgoubi.fai_stepSize0.01cm' u (sqrt($25**2-$29**2)*MeV2eV/c):(($10 - $3)/$3)     w lp pt 50 ps .4 lc 4  tit "  " ,\
    tof(x) axes x1y2 w l lc 2 notit 

     set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_zgoubi.fai_RandTOFvsBrho.eps"  
       replot  
       set terminal X11  
       unset output  

 
pause 1
exit

