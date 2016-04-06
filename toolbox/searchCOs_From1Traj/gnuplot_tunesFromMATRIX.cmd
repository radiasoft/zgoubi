
 set title "Relative time of flight difference"    font "roman,16"   # offset 0,+.7    

 set xlabel "G{/Symbol .g}"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set x2label "E [GeV]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "{/Symbol d}T/T"             font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6

set key top left  font "roman, 12"  samplen 1  

set grid

G =  1.15965218E-03
am = 5.10998920E-01 
c = 2.99792458e8
Emin = 2.5e3
Emax = 15e3
ERef =  7e3  # MeV
pi = 3.14159265359
T_ERef = 9.3300601E-03
GV2MV = 1e-3

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
set x2range  [Emin/1e3:Emax/1e3]
set xrange [G*Emin/am/1e3:G*Emax/am/1e3]
#set yrange [12.81492:12.81515]

 plot \
           "tunesFromMatrix.out" u (G * $12/am):(($10-T_ERef)/T_ERef) axes x1y1 w lp ps .2 lt 8 tit "{/Symbol d}T/T"  ,\
           "tunesFromMatrix.out" u ($12 *GV2MV ):(($10-T_ERef)/T_ERef) axes x2y1 w lp ps .2 lt 8 notit 
 
 set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12
 set output "gnuplot_tunesFromMATRIX_TOF.eps"
 replot
 set terminal X11
 unset output

#################

 set title "Tunes and orbit (at middle of long drift)"    font "roman,16"   # offset 0,+.7    

 set xlabel "G{/Symbol .g}"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set x2label "E [GeV]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "Q_x, Q_y"             font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "x_{orbit} (mm)"          font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics  font "roman,12" nomirror      #offset 0,-.6

set key maxrow 1
set key t c  font "roman, 12"  samplen 1  

set grid

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
set xrange  [*:*]
set x2range [*:*]
set yrange  [*:*]
set y2range [*:*]

 plot \
      "tunesFromMatrix.out" u (G * $12/am):($3) axes x1y1 w l  lt 1 lc 1 tit "Q_x"  ,\
      "tunesFromMatrix.out" u (G * $12/am):($4) axes x1y1 w l  lt 1 lc 3 tit "Q_y"  ,\
      "tunesFromMatrix.out" u ($12 * GV2MV):($1*1e3) axes x2y2 w l  lt 2 tit "x_{co}"  

 set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12
 set output "gnuplot_tunesFromMATRIX_QxQy.eps"
 replot
 set terminal X11
 unset output

      pause 1

 exit
