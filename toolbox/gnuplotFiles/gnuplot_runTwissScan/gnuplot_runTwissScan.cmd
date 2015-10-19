
 set title "Tunes and chromas from Twiss scan"    font "roman,16"   # offset 0,+.7    

 set xlabel "a{/Symbol .g}"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set x2label "E [GeV]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "Q_x, Q_y"             font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "{/Symbol x}_x, {/Symbol x}_y"          font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics  font "roman,12" nomirror      #offset 0,-.6

set key bot left  font "roman, 12"  samplen 1  

set grid

G =  1.15965218E-03
am = 5.10998920E-01 
c = 2.99792458e8
pi = 3.1415926
broRef = 55.038075681301081

#fit_limit=1e-20
#q0 = 0.25
#q1 =1.
#q2 =1.
#q3 =1.
# Qy(x) = q0 + q1/x + q2/x**2 + q3/x**3
# fit [] Qy(x)             "tunesFromMatrix.out" u (($12+am)/am*G):($4) via q0, q1, q2, q3
#print "Qy(44.8339) : ", Qy(44.8339)*138*6
#exit
 
Ggmi = 14 ; Ggma = 49
set xrange  [Ggmi:Ggma] 
set x2range [Ggmi/G*am/1e3:Ggma/G*am/1e3] 
set yrange []
set y2range []

 plot \
      "run_twissScan.out" u (G * $2*broRef*c*1e-6/am):($4) axes x1y1 w lp ps .2 lt 1 lc 1 tit "Q_x"  ,\
      "run_twissScan.out" u (G * $2*broRef*c*1e-6/am):($6) axes x1y1 w lp ps .2 lt 1 lc 3 tit "Q_y"  ,\
      "run_twissScan.out" u (G * $2*broRef*c*1e-6/am):($8) axes x1y2 w lp ps .2 lt 2 lc 1 tit "{/Symbol x}_x"  ,\
      "run_twissScan.out" u (G * $2*broRef*c*1e-6/am):($10) axes x1y2 w lp ps .2 lt 2 lc 3 tit "{/Symbol x}_y"  ,\
      "searchCO.out_COs_12traj"  u (G * $6*broRef*c*1e-6/am):($3+.4) w impulse lt 1 lw .2 lc 2 axes x1y1 tit "Design E" 

#set samples 100000
 set terminal postscript eps blacktext color enh size 8.3cm,6cm "Times-Roman" 12
 set output "gnuplot_runTwissScan_QDQ.eps"
 replot
 set terminal X11
 unset output

      pause 1

 exit
