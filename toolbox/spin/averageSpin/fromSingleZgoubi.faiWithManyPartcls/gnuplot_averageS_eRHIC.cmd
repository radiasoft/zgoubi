
 set title "SPIN DIFFUSION \n Initial zero 6-D emittance"    font "roman,14"   # offset 0,+.7    

 set xlabel "G{/Symbol .g}"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set x2label "E (GeV)"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "Polar. (<cos{/Symbol D f}>)"  font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "{/Symbol s}_{/Symbol f} (deg) \n {/Symbol s}_{E} (MeV) \n G{/Symbol g a} (deg)"  font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics font "roman,12" nomirror      #offset 0,-.6

 set bmargin +4

#set key bottom left  Left font "roman,10"  samplen 1 
set key at graph .53, graph .85  font "roman, 11"  samplen 2 spacin 1.  
 
set grid

G =  1.15965218E-03
am = 5.10998920E-01 
c = 2.99792458e8
p0 =  57.36635309*c
pi = 3.1415926
Einj = 7944


sigSFromzpop = 14.18
sigSFromAvSout = 8.4818
norm =  sigSFromzpop / sigSFromAvSout

fit_limit=1e-20

set xrange  [17:49]
set x2range [17/G*am/1000.:49/G*am/1000.]

 plot \
   "averageS.out_save" u ($1 ):($12 < 12 ? $6       : 1/0)  axes x1y1 w lp pt 6 lc 1 ps .7 tit "Polar"  ,\
   "gnuplot_spinDiff.data" u ($1     ):($4)  axes x1y2 w lp pt 7 lc 3 ps .7 tit "{/Symbol s}_{/Symbol f} (zpop)"  ,\
   "averageS.out_save" u ($1 ):($12 < 12 ? $13*norm : 1/0)  axes x1y2 w lp pt 9 lc 5 ps .7 tit "{/Symbol s}_{/Symbol f} (avrgS)",\
   "theorDiff.out" u ($1 *1e3*G/am):($2)  axes x1y2 w lp pt 11 lc 7 ps .7 tit "{/Symbol s}_{/Symbol f} (th./ring)"  ,\
   "averageS.out_save" u ($1 ):($12 < 12 ? $16 : 1/0)  axes x1y2 w lp pt 9 lc 5 ps .7 tit "{/Symbol s}_{E} (avrgS)" ,\
   "theorDiff.out" u ($1 *1e3*G/am):($3)  axes x1y2 w lp pt 12 lc 10 ps .7 tit "{/Symbol s}_E (th./ring)"  ,\
   "averageS.out_save" u ($1 ):($12 < 12 ? $8/pi*180. : 1/0)  axes x1y2 w lp pt 8 lc 4 ps .7 tit "G{/Symbol ga }  " 


set samples 100000
 set terminal postscript eps blacktext color enh size 7.3cm,5cm "Times-Roman" 12
 set output "gnuplot_averageS_exeydpZro_arcRecent.eps"
 replot
 set terminal X11
 unset output

pause 2
exit
