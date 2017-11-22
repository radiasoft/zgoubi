
 set title "Running average, eRHIC storage ring"    font "roman,16"   # offset 0,+.7    

 set xlabel "s [m]"        font "roman,12"   # offset +4,-.5 rotate by +20  
 set ylabel "{/Symbol bg e}_x [{/Symbol p}m]"             font "roman,12"   #offset -0,-1 rotate by -20 
 set y2label "<S_x>"          font "roman,12"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,10" nomirror
 set x2tics  font "roman,10" nomirror
 set ytics  font "roman,10" nomirror      #offset 0,-.6
 set y2tics  font "roman,10" nomirror      #offset 0,-.6

#set key above maxrows 1 
#set key tm c maxrows 1 width 4
set key t c maxrows 1 width 4
set key font "roman, 12"  samplen 1  

set grid

set samples 10000
#set xrange  []
#set x2range []
#set yrange []

#set y2range [.98:1.]

betx = 0.63309
bety = 0.07606
cm2m=0.01
gma = 18e3/.511

 plot \
      "Run0/runningAverage_fromFai.out" u ($2):(gma*($8 - $7*$7)*cm2m**2/betx) axes x1y1 w l lt 1 linecolor rgb "red" tit "{/Symbol e}_x/{/Symbol p}"  ,\
      "Run1/runningAverage_fromFai.out" u ($2):(gma*($8 - $7*$7)*cm2m**2/betx) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      "Run2/runningAverage_fromFai.out" u ($2):(gma*($8 - $7*$7)*cm2m**2/betx) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      "Run3/runningAverage_fromFai.out" u ($2):(gma*($8 - $7*$7)*cm2m**2/betx) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      "Run0/runningAverage_fromFai.out" u ($2):($19) axes x1y2 w l lt 2 linecolor rgb "blue" tit "<S_x>"  ,\
      "Run1/runningAverage_fromFai.out" u ($2):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit  ,\
      "Run2/runningAverage_fromFai.out" u ($2):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit  ,\
      "Run3/runningAverage_fromFai.out" u ($2):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit    

 set samples 10000
 set terminal postscript eps blacktext color enh size 8.3cm,5cm "Times-Roman" 12
 set output "gnuplot_runAv.eps"
 replot
 set terminal X11
 unset output

pause 2
