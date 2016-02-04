reset 
 
set title "Concentration ellipse surface \n along splitter line \n - Computed from zgoubi.fai -" font "roman,14"

set xlabel "s (m)" font "roman,18"
set ylabel "{/Symbol e}_x, {/Symbol e}_y norm. ({/Symbol m}m)" font "roman,18"
#set y2label "{/Symbol e}_y (m)" font "roman,18"
 
set xtics font "roman,11"
set ytics mirror font "roman,11"
set y2tics nomirror font "roman,11"
 
set key  maxrow 1
set key  top center 
set key font "roman,10" 

#set logscale y 10 
#set logscale y2 10 

am = .511e6
E = 20e9
gam = E/am


#set yrange   [] 
#set y2range  [] 

plot  [] \
    "./lipsFitFromFai_iterate.Out"    u 30:($3*gam*1e6) axes x1y1 w lp ps .4 lt 3 lc 1 tit '{/Symbol e}_x' ,\
    "./lipsFitFromFai_iterate.Out"    u 30:($4*gam*1e6) axes x1y1 w lp ps .4 lt 3 lc 3 tit '{/Symbol e}_y' 

set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12 
 set output "gnuplot_lipsFromFai.eps" 
 replot 
 set terminal X11 
 unset output 


pause 8
exit
##########################################3
##########################################3
##########################################3

