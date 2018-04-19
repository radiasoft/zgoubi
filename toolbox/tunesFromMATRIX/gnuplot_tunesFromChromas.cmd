
# set title "Tunes and chromas from MATRIX scan"    font "roman,16"   # offset 0,+.7    

 set xlabel "m{/Symbol g} [MeV]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "Q_x, Q_y"             font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "{/Symbol x}_x, {/Symbol x}_y"          font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics  font "roman,12" nomirror      #offset 0,-.6

set key t c  font "roman, 12"  samplen 1  
set key maxrow 2


G =  1.15965218E-03
am = 5.10998920E-01 
c = 2.99792458e8
pi = 3.1415926
broRef = 0.14008655052543029

fit_limit=1e-20
q0              = 1.92688      
q1              = -480.79      
q2              = 36102.2      
q3              = -885028      
q4              = 0.99945
q5 = -1
xx = 100.
#q5              = 0.999971    
 xix(x) = q0 + q1*(x-xx) + q2*(x-xx)**2 + q3*(x-xx)**3 + q4*(x-xx)**4  #+ q5*(x-xx)**5
 fit [40:140] xix(x)    "run_twissScan.out" u (($2+1.)*broRef*c*1e-6 < 160.1 ? .9*($2+1.)*broRef*c*1e-6 : 1/0):($10 >-0.5 && $10 <0.3 ? $10 : 1/0 ) via xx, q0, q1, q2, q3, q4  #, q5
print "xi_x(44.8339) : ", xix(44.8339)

r0              = 2.20934      
r1              = -577.69      
r2              = 46918.7      
r3              = -1.27686e+06 
r4              = 0.999884
r5 = -1
yy = 100.
r5              = 0.999997     
 xiy(x) = r0 + r1*(x-yy) + r2*(x-yy)**2 + r3*(x-yy)**3 + r4*(x-yy)**4  + r5*(x-yy)**5
 fit [:150] xiy(x)    "run_twissScan.out" u (($2+1.)*broRef*c*1e-6 < 160.1 ? ($2+1.)*broRef*c*1e-6 : 1/0):($10 >-0.5 && $10 <0.3 ? $10 : 1/0 ) via yy, r0, r1, r2, r3, r4, r5
print "xi_x(44.8339) : ", xix(44.8339)
#exit
 
set xrange  [38:170] 
#set x2range [Ggmi/G*am/1e3:Ggma/G*am/1e3] 
#set yrange []
set y2range [-.8:.3]

set grid

 plot  \
      "run_twissScan.out" u (($2+1.)*broRef*c*1e-6):($4)  axes x1y1 w p pt 5 ps .2  lt 1 lc rgb "red" tit "Q_x"  ,\
      "../../mapModel/lattice/run_twissScan.out" u (($2+1.)*broRef*c*1e-6):($4)  axes x1y1 w p pt 6 ps .2  lt 1 lc rgb "red" notit  ,\
      "run_twissScan.out" u (($2+1.)*broRef*c*1e-6):($7)  axes x1y1 w p pt 5 ps .2  lt 1 lc rgb "blue" tit "Q_y"  ,\
      "../../mapModel/lattice/run_twissScan.out" u (($2+1.)*broRef*c*1e-6):($7)  axes x1y1 w p pt 6 ps .2  lt 1 lc rgb "blue" notit  ,\
      "chromaFromMatrix.out" u (($13)*broRef*c*1e-6):($14*$15 !=0 && $14 <0 && $14 >-.6 ? $14 :1/0)  axes x1y2 w p pt 8 ps .2  lt 1 lc rgb "red" tit "{/Symbol x}_x"  ,\
      "chromaFromMatrix.out" u (($13)*broRef*c*1e-6):($14*$15 !=0 && $15 <0 ? $15 :1/0)  axes x1y2 w p pt 8 ps .2  lt 1 lc rgb "blue" tit "{/Symbol x}_y"   ,\
       "orbits_4.data"  u ($6*broRef*c*1e-6):($3+.4) w impulse lt 1 lw 2 lc rgb "green" axes x1y1 tit "Design E" 
       
# plot  \
#      "run_twissScan.out" u (($2+1.)*broRef*c*1e-6):($4)  axes x1y1 smooth csplines w l lw 2 lt 1 lc rgb "red" tit "Q_x"  ,\
#      "../../mapModel/lattice/run_twissScan.out" u (($2+1.)*broRef*c*1e-6):($4)  axes x1y1 smooth csplines w l lw 2 lt 1 lc rgb "red" notit  ,\
#      "run_twissScan.out" u (($2+1.)*broRef*c*1e-6):($7)  axes x1y1 smooth csplines w l lw 2 lt 1 lc rgb "blue" tit "Q_y"  ,\
#      "../../mapModel/lattice/run_twissScan.out" u (($2+1.)*broRef*c*1e-6):($7)  axes x1y1 smooth csplines w l lw 2 lt 1 lc rgb "blue" notit  ,\
#        (x>=40 && x<=152 ? xix(x) :1/0) axes x1y2 lw 1 lc rgb "red" tit "{/Symbol x}_x" ,\
#        (x>=40 && x<=152 ? xiy(x) :1/0) axes x1y2 lw 1 lc rgb "blue" tit "{/Symbol x}_y"  ,\
#       "orbits_4.data"  u ($6*broRef*c*1e-6):($3+.4) w impulse lt 1 lw 2 lc rgb "green" axes x1y1 tit "Design E" 
       
#set samples 100000
 set terminal postscript eps blacktext color enh size 8.3cm,6cm "Times-Roman" 12
 set output "gnuplot_QsAndChromas.eps"
 replot
 set terminal X11
 unset output

      pause 1

 exit
