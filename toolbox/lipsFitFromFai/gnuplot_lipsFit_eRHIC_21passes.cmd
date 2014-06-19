reset

set title "EMITTANCES" font "roman,14"
 
set xlabel "G{/Symbol g}" font "roman,13"
set x2label "E (GeV)" font "roman,13" 
#set ylabel "{/Symbol e}_x, {/Symbol e}_y, geom. (nm)" font "roman,12"
set ylabel  "{/Symbol e}_x, {/Symbol e}_y, norm. ({/Symbol m}m)" font "roman,12"
#set y2label "{/Symbol e}_x, {/Symbol e}_y, norm. ({/Symbol m}m)" font "roman,12"
 
set xtics nomirror font "roman,11"
set x2tics nomirror font "roman,11" 
set ytics mirror font "roman,11"
set y2tics nomirror font "roman,11"
 
set key font "roman,10"  top center samplen 2 spacing .9  

set grid 

#set logscale y 10 
#set logscale y2 10 

am = 0.51099892
G = 1.15965218076e-3

set xrange  [17:49]
set x2range [17/G*am/1000.:49/G*am/1000.]
#set yrange [-0.5:50]
#set y2range [-0.05:.7]

plot \
           "orbits_11.data" u ($8/am*G):($3+60) w impulse lt 1 lw .8 lc 3 axes x1y1 tit "     design Energy" ,\
    "./lipsFitFromFai_iterate_save.Out"    u (G*($13/am)):($3*1e6 * ($13/am)) axes x1y1 w lp ps .8  pt 6 lc 1 tit '{/Symbol e}_{x,norm.}' ,\
    "./lipsFitFromFai_iterate_save.Out"    u (G*($13/am)):($4*1e6 * ($13/am)) axes x1y1 w lp ps .8  pt 8 lc 3 tit '{/Symbol e}_{y,norm.}'


#plot \
#           "orbits_11.data" u ($8/am*G):($3+4) w impulse lt 1 lw 1.2 lc 3 axes x1y1 tit "     design Energy" ,\
#    "./lipsFitFromFai_iterate_save.Out"    u ($13/1e3)   :($3*1e9)            axes x2y1 w lp ps .8  pt 5 lc 1 tit '{/Symbol e}_x, geom' ,\
#    "./lipsFitFromFai_iterate_save.Out"    u (G*($13/am)):($3*1e6 * ($13/am)) axes x1y2 w lp ps .8  pt 6 lc 1 tit '{/Symbol e}_{x,norm.}' ,\
#    "./lipsFitFromFai_iterate_save.Out"    u ($13/1e3)   :($4*1e9)            axes x2y1 w lp ps .8  pt 7 lc 3 tit '{/Symbol e}_y, geom' ,\
#    "./lipsFitFromFai_iterate_save.Out"    u (G*($13/am)):($4*1e6 * ($13/am)) axes x1y2 w lp ps .8  pt 8 lc 3 tit '{/Symbol e}_{y,norm.}'


set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_emittances_linac1.322_CAVITESRComp_21passes.eps" 
 replot 
 set terminal X11 
 unset output 

pause 2
exit
##########################################3
##########################################3
##########################################3

