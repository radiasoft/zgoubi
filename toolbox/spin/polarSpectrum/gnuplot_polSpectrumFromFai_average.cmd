
 set title "Polarization vs. time at store, 18 GeV \n From polSpectrumFromFai_average.res"    font "roman,16"   # offset 0,+.7    

 set xlabel "{/Symbol t} [{/Symbol }10^3 turn]"        font "roman,12"   # offset +4,-.5 rotate by +20  
 set x2label "{/Symbol t} [{/Symbol } turn])"        font "roman,12"   # offset +4,-.5 rotate by +20  
 set ylabel "<S_l>"          font "roman,12"   #offset -0,-1 rotate by -20 
 set y2label "<S_l>"          font "roman,12"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,10" nomirror
 set x2tics  font "roman,10" nomirror
 set ytics  font "roman,10" nomirror      #offset 0,-.6
 set y2tics  font "roman,10" nomirror      #offset 0,-.6

set key t r maxcol 1 # width 1
set key font "roman, 10"  samplen 1  

set grid
set samples 10000 

 gma = 18e3/.511
 frev = 3e8/3835.

     P0 = 1.
     tauD = 1e6;  tauD1 = 1./tauD
     P(x) = P0*exp(-x * tauD1)
     FIT_LIMIT = 1e-6
     fit P(x) "polSpectrumFromFai_average.res" u ($4):($5) via tauD1
     tauD = 1./tauD1

tauST = 30.*60 * frev  # turn
tauST1 = 1./ tauST
tauEq = 1./ (tauD1 + tauST1)

     print ' tau_diffusion  (turn) = ',tauD,' turn,  ',tauD/frev/60.,' minutes'
     print ' tau_eq  (turn) = ',tauEq,' turn,  ',tauEq/frev/60.,' minutes'

set x2range [0:2*tauD]
tauRS=500

 plot \
     "polSpectrumFromFai_average.res" u ($4 < 1000*tauRS ? $4/1e3 : 1/0):($5) axes x1y1 w lp ps .3 pt 5  lt rgbcolor  "red" tit "x1y1: S_l zoom" ,\
     "polSpectrumFromFai_average.res" u ($4):($5) axes x2y2 w lp ps .3 pt 4  lt rgbcolor  "blue" tit "x2y2: S_l" ,\
     P(x) axes x2y2 tit "x2y2: e^{-t/tauEq}"

 set samples 10000
 set terminal postscript eps blacktext color enh "Times-Roman" 12
 set output "gnuplot_polSpectrumFromFai_average.eps"
 replot
 set terminal X11
 unset output
  
pause 2
exit

  

