set terminal postscript eps blacktext color enh size 15cm,7cm "Times-Roman" 17
set output "twiss.eps"

set title"GENERALIZED TWISS' PARAMETERS"
set xlabel"Arclength (in m)"
#set xrange[0:200]
set ylabel"{/Symbol b}_I , {/Symbol b}_{II} (in m)"
#set y2label"{/Symbol a}_I , {/Symbol a}_{II}"

set ytics nomirror
set y2tics nomirror
set grid
set key outside

plot [][] \
 "twiss.out" u ($3):($4) w lp pt 1 ps 1. tit "{/Symbol b}_I", \
 "twiss.out" u ($3):($5) w lp pt 2 ps 1. tit "{/Symbol b}_{II}"#, \
# "twiss.out" u ($1):($4) axes x1y2 w lp pt 1 ps 1. tit "{/Symbol a}_I", \
# "twiss.out" u ($1):($5) axes x1y2 w lp pt 2 ps 1. tit "{/Symbol a}_{II}"# 


####################################################################################


set terminal postscript eps blacktext color enh size 15cm,7cm "Times-Roman" 17
set output "rPARAM.eps"

set title""
set title"COUPLING PARAMETER r"
set xlabel"Arclength (in m)"
#set xrange[0:200]
set ylabel"r"
#set yrange[0.7:1.1]
set y2label""
#set y2range[0.7:1.1]

set ytics mirror
set grid
set key inside
#set logscale y  
#set logscale y2 

plot [][] \
 "twiss.out" u ($3):($10) w lp pt 1 ps 1. tit "r"#, \


exit
