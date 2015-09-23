
 set title "SPIN DIFFUSION \n Initial Gaussian {/Symbol e}_x{/Symbol \273 e}_y{/Symbol \273}50{/Symbol p} rms, {/Symbol s}_{dp/p}=0"    font "roman,14"   # offset 0,+.7    

 set xlabel "a{/Symbol .g}"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set x2label "E (GeV)"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "Polar. (<cos{/Symbol D f}>)"  font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "{/Symbol s}_{/Symbol f} (deg.) \n a{/Symbol g a} (deg.)"             font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics font "roman,12" nomirror      #offset 0,-.6

 set bmargin +4

#set key bottom left  Left font "roman,10"  samplen 1 
set key center left  font "roman, 13"  samplen 2 spacin .7  
 
set grid

a =  1.15965218E-03
am = 5.10998920E-01 
c = 2.99792458e8
p0 =  57.36635309*c
pi = 3.1415926
Einj = 7944


sigSFromzpop = 1.464304E+01           # var=28,26
sigSFromAvSout = 8.71018952E+00
norm = sigSFromzpop / sigSFromAvSout #sigSFromzpop / sigSFromAvSout

fit_limit=1e-20

set xrange  [17:49]
set x2range [17/a*am/1000.:49/a*am/1000.]

 plot \
<<<<<<< .mine
   "averageS.out_save" u ($1     ):($12 < 12 ? $6       : 1/0)  axes x1y1 w lp pt 6 lc 1 ps .7 tit "<Polar>"  ,\
   "averageS.out_save" u ($1     ):($12 < 12 ? $13*norm : 1/0)  axes x1y2 w lp pt 7 lc 3 ps .7 tit "{/Symbol s}_{/Symbol f}" ,\
   "averageS.out_save" u ($1     ):($12 < 12 ? $8/pi*180. : 1/0)  axes x1y2 w lp pt 8 lc 4 ps .7 tit "a{/Symbol ga }" 
=======
   "averageS.out_save" u ($1     ):($12 < 12 ? $6       : 1/0)  axes x1y1 w lp pt 6 lc 1 ps .7 tit "<Polar>"  ,\
   "averageS.out_save" u ($1     ):($12 < 12 ? $13*norm : 1/0)  axes x1y2 w lp pt 7 lc 3 ps .7 tit "{/Symbol s}_{/Symbol f}" ,\
   "averageS.out_save" u ($1     ):($12 < 12 ? $8/pi*180. : 1/0)  axes x1y2 w lp pt 8 lc 4 ps .7 tit "G{/Symbol ga }" 
>>>>>>> .r572


set samples 100000
 set terminal postscript eps blacktext color enh size 7.3cm,5cm "Times-Roman" 12
 set output "gnuplot_averageS_exey50-dp0_arcRecent.eps"
 replot
 set terminal X11
 unset output

pause 2
exit
