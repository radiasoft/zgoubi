set title "Running average, eRHIC storage ring"    font "roman,16"   # offset 0,+.7    

 set xlabel "s [m]"        font "roman,12"   # offset +4,-.5 rotate by +20  
 set ylabel "{/Symbol bg e}_x [{/Symbol p m}m]"             font "roman,12"   #offset -0,-1 rotate by -20 
 set y2label "<S_x>"          font "roman,12"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,10" nomirror
 set x2tics  font "roman,10" nomirror
 set ytics  font "roman,10" nomirror      #offset 0,-.6
 set y2tics  font "roman,10" nomirror      #offset 0,-.6

set logscale y 10

#set key above maxrows 1 
#set key tm c maxrows 1 width 4
set key b c maxrows 1 width 4
set key font "roman, 12"  samplen 1  

set grid
set samples 10000 

betx = 0.63309
bety = 0.07606
cm2m=0.01
mrd2rd=0.001
m2mum = 1e6
gma = 18e3/.511

set xrange  [1100:]
#set x2range [:]
set yrange [:]
set y2range [:1.05]

 plot \
      "Run000/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($7*cm2m)**2/betx + betx*($9*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" tit "{/Symbol e}_x"  ,\
      for [i=1:9] "Run00".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($7*cm2m)**2/betx + betx*($9*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run01".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($7*cm2m)**2/betx + betx*($9*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run02".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($7*cm2m)**2/betx + betx*($9*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run03".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($7*cm2m)**2/betx + betx*($9*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run04".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($7*cm2m)**2/betx + betx*($9*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run05".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($7*cm2m)**2/betx + betx*($9*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run06".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($7*cm2m)**2/betx + betx*($9*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run07".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($7*cm2m)**2/betx + betx*($9*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run08".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($7*cm2m)**2/betx + betx*($9*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run09".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($7*cm2m)**2/betx + betx*($9*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run10".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($7*cm2m)**2/betx + betx*($9*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      for [i=0:9] "Run11".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($7*cm2m)**2/betx + betx*($9*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" notit  ,\
      "Run000/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" tit "<S_x>"  ,\
      for [i=1:9] "Run00".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run01".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run02".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run03".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run04".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run05".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run06".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run07".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run08".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run09".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run10".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run11".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" notit 

 set samples 10000
 set terminal postscript eps blacktext color enh size 8.3cm,5cm "Times-Roman" 12
 set output "gnuplot_runAvEpsx.eps"
 replot
 set terminal X11
 unset output

pause 2
set ylabel "{/Symbol bg e}_y [{/Symbol p m}m]"             font "roman,12"   #offset -0,-1 rotate by -20 

 plot \
      "Run000/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($11*cm2m)**2/bety + bety*($13*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" tit "{/Symbol e}_y"  ,\
      for [i=1:9] "Run00".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($11*cm2m)**2/bety + bety*($13*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" tit "{/Symbol e}_y"  ,\
      for [i=0:9] "Run01".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($11*cm2m)**2/bety + bety*($13*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" tit "{/Symbol e}_y"  ,\
      for [i=0:9] "Run02".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($11*cm2m)**2/bety + bety*($13*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" tit "{/Symbol e}_y"  ,\
      for [i=0:9] "Run03".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($11*cm2m)**2/bety + bety*($13*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" tit "{/Symbol e}_y"  ,\
      for [i=0:9] "Run04".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($11*cm2m)**2/bety + bety*($13*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" tit "{/Symbol e}_y"  ,\
      for [i=0:9] "Run05".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(m2mum*gma*(($11*cm2m)**2/bety + bety*($13*mrd2rd)**2)) axes x1y1 w l lt 1 linecolor rgb "red" tit "{/Symbol e}_y"  ,\
      "Run000/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" tit "<S_x>"  ,\
      for [i=1:9] "Run00".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" tit "<S_x>"  ,\
      for [i=0:9] "Run01".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" tit "<S_x>"  ,\
      for [i=0:9] "Run02".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" tit "<S_x>"  ,\
      for [i=0:9] "Run03".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" tit "<S_x>"  ,\
      for [i=0:9] "Run04".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" tit "<S_x>"  ,\
      for [i=0:9] "Run05".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($19) axes x1y2 w l lt 2 linecolor rgb "blue" tit "<S_x>"  

 set samples 10000
 set terminal postscript eps blacktext color enh size 8.3cm,5cm "Times-Roman" 12
 set output "gnuplot_runAvEpsy.eps"
 replot
 set terminal X11
 unset output


pause 2


unset logscale y 

set ylabel "<x> [m]"             font "roman,12"   #offset -0,-1 rotate by -20 
set y2label "<y> [m]"             font "roman,12"   #offset -0,-1 rotate by -20 

 plot \
      "Run000/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(($7*cm2m)) axes x1y1 w l lt 1 linecolor rgb "red" tit "x"  ,\
      for [i=1:9] "Run00".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(($7*cm2m)) axes x1y1 w l lt 1 linecolor rgb "red" notit ,\
      for [i=0:9] "Run01".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(($7*cm2m)) axes x1y1 w l lt 1 linecolor rgb "red" notit ,\
      for [i=0:9] "Run02".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(($7*cm2m)) axes x1y1 w l lt 1 linecolor rgb "red" notit ,\
      for [i=0:9] "Run03".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(($7*cm2m)) axes x1y1 w l lt 1 linecolor rgb "red" notit ,\
      for [i=0:9] "Run04".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(($7*cm2m)) axes x1y1 w l lt 1 linecolor rgb "red" notit ,\
      for [i=0:9] "Run05".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):(($7*cm2m)) axes x1y1 w l lt 1 linecolor rgb "red" notit ,\
      "Run000/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($11) axes x1y2 w l lt 2 linecolor rgb "blue" tit "y"  ,\
      for [i=1:9] "Run00".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($11) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run01".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($11) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run02".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($11) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run03".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($11) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run04".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($11) axes x1y2 w l lt 2 linecolor rgb "blue" notit ,\
      for [i=0:9] "Run05".i."/runningAverage_fromFai.out" u (int($2/100)*100==$2 ? $2 : 1/0):($11) axes x1y2 w l lt 2 linecolor rgb "blue" notit 

 set samples 10000
 set terminal postscript eps blacktext color enh size 8.3cm,5cm "Times-Roman" 12
 set output "gnuplot_runAvxy.eps"
 replot
 set terminal X11
 unset output

pause 88 
