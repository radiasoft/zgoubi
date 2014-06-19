
 set title "SPIN DIFFUSION, {/Symbol e}_x=50{/Symbol p} (SR loss compensation by linac)"    font "roman,12"   # offset 0,+.7    

 set xlabel "Design G{/Symbol .g}"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set x2label "Design energy [GeV]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "{/Symbol s}_{/Symbol f} [deg.]"             font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "Polar. (<cos{/Symbol D f}>)"  font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics font "roman,12" nomirror      #offset 0,-.6

#set key bottom left  Left font "roman,10"  samplen 1 
set key center left  font "roman, 12"  samplen 2  
 
set grid

G =  1.15965218E-03
am = 5.10998920E-01 
c = 2.99792458e8
p0 =  57.36635309*c
pi = 3.1415926
Einj = 7944

norm = 11.8393845 / 7.30514979

fit_limit=1e-20

set xrange  [17:49]
set x2range [17/G*am/1000.:49/G*am/1000.]

 plot \
     "averageS.out" u ($1     ):($12 < 12 ?  $6 : 1/0)  axes x1y2 w lp ps .8 pt 10 lc 3 notit ,\
     "averageS.out" u ($1     ):($12 < 12 ?  $6 : 1/0)  axes x1y2 w lp ps 1. pt 12 lc 3 tit "Polar. (<cos{/Symbol D f}>)" ,\
     "averageS.out" u ($1/G*am):($12 < 12 ? $13*norm : 1/0)  axes x2y1 w lp ps .8 pt 6 notit ,\
     "averageS.out" u ($1     ):($12 < 12 ? $13*norm : 1/0)  axes x1y1 w lp ps 1. pt 8 lc 1 tit "{/Symbol s}_{/Symbol f}"  


set samples 100000
 set terminal postscript eps blacktext color enh size 7.3cm,5cm "Times-Roman" 12
 set output "gnuplot_averageS_exey50_ZroDfct.eps"
 replot
 set terminal X11
 unset output

pause 222
exit
