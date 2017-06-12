
set key maxcol 1
set key t l

#set logscale y 

set tit 'TOF vs. particle rigidity, from zgoubi.fai'

set xtics nomirror
set x2tics nomirror
set ytics 

set xlabel 'B{/Symbol r}   /T.m'
set x2label 'kin. E   /MeV'
set ylabel 'time of flight   /{/Symbol m}s'

# electron : 
#am = 0.511
# proton :
am = 938.27203
m2cm = 100.
MeV2eV = 1e6
c = 2.99792458e8
B = 0.5   # [T]
rho(x) = m2cm * x/B    # [cm]
pi = 4. * atan(1.)
Nsector = 6
s2mus = 1e6

# Plotting thoretical TOF establishes that slight TOF increase observed from tracking is as expected from relativity. 
TOF(x) = (2.*pi* x/B /Nsector) / ((x*c)/sqrt((x*c)**2 + (am*MeV2eV)**2)) /c*s2mus  # t = Circ / v

Ekmi = 0.2
Ekma = 5
Emi = Ekmi + am
Ema = Ekma + am
Brmi = sqrt(Emi**2-am**2)*MeV2eV/c
Brma = sqrt(Ema**2-am**2)*MeV2eV/c

set xrange [Brmi:Brma]
set x2range [Ekmi:Ekma]

plot  \
   'zgoubi.fai_stepSize0.1cm' u (sqrt($25**2-$29**2)*MeV2eV/c):15 w lp pt 5 ps .4 lc 1  tit "step 0.1cm" ,\
   'zgoubi.fai_stepSize0.03cm' u (sqrt($25**2-$29**2)*MeV2eV/c):15 w lp pt 5 ps .4 lc 1  tit "step 0.03cm" ,\
   'zgoubi.fai_stepSize0.01cm' u (sqrt($25**2-$29**2)*MeV2eV/c):15 w lp pt 5 ps .4 lc 3  tit "step 0.01cm"  ,\
    TOF(x) w l lc 4 tit 'R(B{/Symbol r}), theor.' 
     set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_zgoubi.fai_TOFvsBrho_3steps.eps"  
       replot  
       set terminal X11  
       unset output  
 
pause 1


set xlabel 'R   /cm'
set y2label 'Orbit angle   /rad'

set ytics nomirror
set y2tics nomirror

set xrange [Brmi/B*m2cm:Brma/B*m2cm]
set x2range [Ekmi:Ekma]


plot  \
   'zgoubi.fai_stepSize0.1cm' u 10:15 w lp pt 5 ps .4 lc 1  tit "step 0.1cm" ,\
   'zgoubi.fai_stepSize0.03cm' u 10:15 w lp pt 5 ps .4 lc 1  tit "step 0.03cm" ,\
   'zgoubi.fai_stepSize0.01cm' u 10:15 w lp pt 5 ps .4 lc 3  tit "step 0.01cm" ,\
   'zgoubi.fai_stepSize0.01cm' u 10:(abs($11) < .0005 ? $11 : 1/0) axes x1y2 w lp  ps .4  tit "angle"

#plot  \
#    TOF(x) w l lc 1 tit 'R(B{/Symbol r}), theor.' 

     set terminal postscript eps blacktext color  enh size 8.3cm,4cm "Times-Roman" 12  
       set output "gnuplot_zgoubi.fai_TOFvsR.eps"  
       replot  
       set terminal X11  
       unset output  

 
pause 1
exit

