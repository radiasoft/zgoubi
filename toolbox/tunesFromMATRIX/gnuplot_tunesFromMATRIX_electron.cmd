
# set title "Time of flight"    font "roman,16"   # offset 0,+.7    

 set xlabel "a{/Symbol .g}"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set x2label "E [GeV]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "dT/T_{ref}, dL/L_{ref}"             font "roman,16"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics  font "roman,12" nomirror      #offset 0,-.6

set key top left  font "roman, 12"  samplen 1  

set grid

G =  1.15965218E-03
am = 5.10998920E-01 
c = 2.99792458e8
pi = 3.1415926

#fit_limit=1e-20
#q0 = 0.25
#q1 =1.
#q2 =1.
#q3 =1.
# Qy(x) = q0 + q1/x + q2/x**2 + q3/x**3
# fit [] Qy(x)             "tunesFromMatrix.out" u (($12+am)/am*G):($4) via q0, q1, q2, q3
#print "Qy(44.8339) : ", Qy(44.8339)*138*6
#exit

#set samples 100000
#set xrange  [17:49]
#set x2range [17/G*am/1000.:49/G*am/1000.]
###set yrange [12.81492:12.81515]


Ekref = 41.989e-3
Lref = 44.507374  # cm
Tref = 1.48471352E-03 # mu_s

 plot \
           "tunesFromMatrix_31trj.out_save" u (G * $12/am):(($11-Lref)/Lref) w lp ps .2 lt 8 lc 1 tit "dL/L_{ref}"  ,\
           "tunesFromMatrix_31trj.out_save" u (G * $12/am):(($10-Tref)/Tref) w lp ps .2 lt 8 lc 3 tit "dT/T_{ref}"  

 set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12
 set output "gnuplot_tunesFromMATRIX_pathAndTOF.eps"
 replot
 set terminal X11
 unset output

pause 1

#################

 set title "Tunes"    font "roman,16"   # offset 0,+.7    

 set xlabel "a{/Symbol .g}"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set x2label "E [GeV]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "Q_x, Q_y"             font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "x_{orbit} (mm)"          font "roman,13"   #offset -0,-1 rotate by -20 

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

#fit_limit=1e-20
#q0 = 0.25
#q1 =1.
#q2 =1.
#q3 =1.
# Qy(x) = q0 + q1/x + q2/x**2 + q3/x**3
# fit [] Qy(x)             "tunesFromMatrix.out" u (($12+am)/am*G):($4) via q0, q1, q2, q3
#print "Qy(44.8339) : ", Qy(44.8339)*138*6
#exit

#set samples 100000
#set xrange  [17:49]
#set x2range [17/G*am/1000.:49/G*am/1000.]
#set yrange [*:*]
#set y2range [-.5:.5]

 plot \
      "tunesFromMatrix_31trj.out_save" u (G * $12/am):($3) axes x1y1 w lp ps .2 lt 20 tit "Q_x"  ,\
      "tunesFromMatrix_31trj.out_save" u (G * $12/am):($4) axes x1y1 w lp ps .2 lt 30 tit "Q_y"  ,\
      "tunesFromMatrix_31trj.out_save" u (G * $12/am):($1*1e3) axes x1y2 w lp ps .2 lt 40 tit "Orbit"  

#      "tunesFromMatrix_31trj.out_save" u (G * $12/am):(G * $12/am - int(G * $12/am)) axes x1y1 w lp ps .2 lt 30 tit "frac(a.{/Symbol g})"  ,\

set samples 100000
 set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12
 set output "gnuplot_tunesFromMATRIX_QxQy.eps"
 replot
 set terminal X11
 unset output

      pause 1

 exit
