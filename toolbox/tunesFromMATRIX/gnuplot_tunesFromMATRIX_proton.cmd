
# set title "Time of flight"    font "roman,16"   # offset 0,+.7    

 set xlabel "G{/Symbol .g}"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set x2label "E [MeV]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "dT/T_{ref}, dL/L_{ref}"             font "roman,16"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics  font "roman,12" nomirror      #offset 0,-.6

set key top left  font "roman, 12"  samplen 1  

set grid

G =  1.79284735
am = 938.27203 
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
           "FFAG/tunesFromMatrix.out_FFAG" u (G * $12/am):(($10-Tref)/Tref) w lp ps .2 lt 8 lc 1 tit "dT/T_{ref}"   ,\
           "IX-567.IR-41.IZ-5/tunesFromMatrix.out_IX-567.IR-41.IZ-5" u (G * $12/am):(($10-Tref)/Tref) w lp ps .2 lt 8 lc 3 tit "dT/T_{ref}"   ,\
           "IX-567.IR-81.IZ-9/tunesFromMatrix.out_IX-567.IR-81.IZ-9" u (G * $12/am):(($10-Tref)/Tref) w lp ps .2 lt 8 lc 3 tit "dT/T_{ref}"   

 set terminal postscript eps blacktext color enh size 8.3cm,4cm "Times-Roman" 12
 set output "gnuplot_tunesFromMATRIX_pathAndTOF.eps"
 replot
 set terminal X11
 unset output

pause 1

#################

 set title "Tunes"    font "roman,16"   # offset 0,+.7    

 set xlabel "G{/Symbol .g}"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set x2label "E [MeV]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "Q_x, Q_y"             font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "x_{orbit} (cm)"          font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics  font "roman,12" nomirror      #offset 0,-.6

set grid

G =  1.79284735
am = 938.27203 
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
set xrange  [(12+am)/am*G:(151+am)/am*G]
set x2range [12:151]
#set yrange [*:*]
#set y2range [-.5:.5]

set key maxrow 4
set key c r  font "roman, 10"  samplen 1  

  plot \
       "./FFAG/tunesFromMatrix.out_FFAG" u (G * ($12+am)/am):($3>0 ? $3 : 1/0) axes x1y1 w lp ps .2 lt 20 tit "FFAG : Q_x"  ,\
       "./IX-567.IR-41.IZ-5/tunesFromMatrix.out" u (G * ($12+am)/am):($3>0 ? $3 : 1/0) axes x1y1 w lp ps .7 lt 21 tit "Map/2cm : Q_x"  ,\
       "./IX-567.IR-81.IZ-9/tunesFromMatrix.out_IX-567.IR-81.IZ-9" u (G * ($12+am)/am):($3>0 ? $3 : 1/0) axes x1y1 w lp ps .7 lt 22 tit "Map/1cm :Q_x"  ,\
       "./IX-1001.IR-161.IZ-41/tunesFromMatrix.out" u (G * ($12+am)/am):($3>0 ? $3 : 1/0) axes x1y1 w lp ps .7 lt 23 tit "Map/5mm :Q_x"  ,\
       "./FFAG/tunesFromMatrix.out_FFAG" u (G * ($12+am)/am):($4>0 ? $4 : 1/0) axes x1y1 w lp ps .2 lt 30 tit "Q_y"  ,\
       "./IX-567.IR-41.IZ-5/tunesFromMatrix.out" u (G * ($12+am)/am):($4>0 ? $4 : 1/0) axes x1y1 w lp ps .7 lt 31 tit "Q_y"  ,\
       "./IX-567.IR-81.IZ-9/tunesFromMatrix.out_IX-567.IR-81.IZ-9" u (G * ($12+am)/am):($4>0 ? $4 : 1/0) axes x1y1 w lp ps .7 lt 32 tit "Q_y"  ,\
       "./IX-1001.IR-161.IZ-41/tunesFromMatrix.out" u (G * ($12+am)/am):($4>0 ? $4 : 1/0) axes x1y1 w lp ps .7 lt 33 tit "Q_y"  ,\
       "./FFAG/tunesFromMatrix.out_FFAG" u (G * ($12+am)/am):($1*1e2) axes x1y2 w lp ps .2 lt 40 tit "Orbit"   ,\
       "./IX-567.IR-41.IZ-5/tunesFromMatrix.out" u (G * ($12+am)/am):($1*1e2) axes x1y2 w lp ps .7 lt 41 tit "Orbit"  ,\
       "./IX-567.IR-81.IZ-9/tunesFromMatrix.out_IX-567.IR-81.IZ-9" u (G * ($12+am)/am):($1*1e2) axes x1y2 w lp ps .7 lt 42 tit "Orbit"  ,\
       "./IX-1001.IR-161.IZ-41/tunesFromMatrix.out" u (G * ($12+am)/am):($1*1e2) axes x1y2 w lp ps .7 lt 43 tit "Orbit"  

#set key maxrow 2
#set key c r  font "roman, 10"  samplen 1  

# plot \
#      "./FFAG/tunesFromMatrix.out_FFAG" u (G * ($12+am)/am):($3>0 ? $3 : 1/0) axes x1y1 w lp ps .2 lt 20 tit "FFAG : Q_x"  ,\
#      "./IX-1001.IR-161.IZ-41/tunesFromMatrix.out" u (G * ($12+am)/am):($3>0 ? $3 : 1/0) axes x1y1 w lp ps .7 lt 23 tit "Map/5mm :Q_x"  ,\
#      "./FFAG/tunesFromMatrix.out_FFAG" u (G * ($12+am)/am):($4>0 ? $4 : 1/0) axes x1y1 w lp ps .2 lt 30 tit "Q_y"  ,\
#      "./IX-1001.IR-161.IZ-41/tunesFromMatrix.out" u (G * ($12+am)/am):($4>0 ? $4 : 1/0) axes x1y1 w lp ps .7 lt 33 tit "Q_y"  ,\
#      "./FFAG/tunesFromMatrix.out_FFAG" u (G * ($12+am)/am):($1*1e3) axes x1y2 w lp ps .2 lt 40 tit "Orbit"   ,\
#      "./IX-1001.IR-161.IZ-41/tunesFromMatrix.out" u (G * ($12+am)/am):($1*1e3) axes x1y2 w lp ps .7 lt 43 tit "Orbit"  

 set terminal postscript eps blacktext color enh size 8.3cm,7cm "Times-Roman" 12
 set output "gnuplot_tunesFromMATRIX_QxQy.eps"
 replot
 set terminal X11
 unset output

      pause 1

 exit
