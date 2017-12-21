set title "Running average, final polarization state \n eRHIC storage ring"    font "roman,16"   # offset 0,+.7    

 set xlabel "a{/Symbol g}"        font "roman,12"   # offset +4,-.5 rotate by +20  
 set ylabel "(SZmax+SZmin)/2"             font "roman,12"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,10" mirror
 set ytics  font "roman,10" mirror      #offset 0,-.6

#set key above maxrows 1 
#set key tm c maxrows 1 width 4
set key t c maxrows 1 width 4
set key font "roman, 12"  samplen 1  

set grid
set samples 10000 

betx = 0.63309
bety = 0.07606
cm2m=0.01
mrd2rd=0.001
m2mum = 1e6
gma = 18e3/.511

#set xrange  [1100:]
#set x2range [:]

set yrange [:1.05]

lastTurn=1000

 plot \
      "Run0/runningAverage_fromFai.out" u ($2==lastTurn ? $34/$35*$37 : 1/0):(0.5*($31+$32)) axes x1y1 w lp ps .4 lt 1 linecolor rgb "red" tit "(SZmi+SZma)/2"  ,\
      for [i=1:9] "Run".i."/runningAverage_fromFai.out" u ($2==lastTurn ? $34/$35*$37 : 1/0):(0.5*($31+$32)) axes x1y1 w lp ps .4 lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run1".i."/runningAverage_fromFai.out" u ($2==lastTurn ? $34/$35*$37 : 1/0):(0.5*($31+$32)) axes x1y1 w lp ps .4 lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run2".i."/runningAverage_fromFai.out" u ($2==lastTurn ? $34/$35*$37 : 1/0):(0.5*($31+$32)) axes x1y1 w lp ps .4 lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run3".i."/runningAverage_fromFai.out" u ($2==lastTurn ? $34/$35*$37 : 1/0):(0.5*($31+$32)) axes x1y1 w lp ps .4 lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run4".i."/runningAverage_fromFai.out" u ($2==lastTurn ? $34/$35*$37 : 1/0):(0.5*($31+$32)) axes x1y1 w lp ps .4 lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run5".i."/runningAverage_fromFai.out" u ($2==lastTurn ? $34/$35*$37 : 1/0):(0.5*($31+$32)) axes x1y1 w lp ps .4 lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run6".i."/runningAverage_fromFai.out" u ($2==lastTurn ? $34/$35*$37 : 1/0):(0.5*($31+$32)) axes x1y1 w lp ps .4 lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run7".i."/runningAverage_fromFai.out" u ($2==lastTurn ? $34/$35*$37 : 1/0):(0.5*($31+$32)) axes x1y1 w lp ps .4 lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run8".i."/runningAverage_fromFai.out" u ($2==lastTurn ? $34/$35*$37 : 1/0):(0.5*($31+$32)) axes x1y1 w lp ps .4 lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run9".i."/runningAverage_fromFai.out" u ($2==lastTurn ? $34/$35*$37 : 1/0):(0.5*($31+$32)) axes x1y1 w lp ps .4 lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run10".i."/runningAverage_fromFai.out" u ($2==lastTurn ? $34/$35*$37 : 1/0):(0.5*($31+$32)) axes x1y1 w lp ps .4 lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run11".i."/runningAverage_fromFai.out" u ($2==lastTurn ? $34/$35*$37 : 1/0):(0.5*($31+$32)) axes x1y1 w lp ps .4 lt 1 linecolor rgb "red" notit   ,\
      for [i=0:7] "Run12".i."/runningAverage_fromFai.out" u ($2==lastTurn ? $34/$35*$37 : 1/0):(0.5*($31+$32)) axes x1y1 w lp ps .4 lt 1 linecolor rgb "red" notit  

 set samples 10000
 set terminal postscript eps blacktext color enh size 8.3cm,5cm "Times-Roman" 12
 set output "gnuplot_runAv_polFinal.eps"
 replot
 set terminal X11
 unset output

exit
