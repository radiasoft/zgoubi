#plot "my.dat" every A:B:C:D:E:F
#A: line increment
#B: data block increment
#C: The first line        (line count starts at 0 !!)
#D: The first data block
#E: The last line
#F: The last data block

 set title "Spin rotation by RHIC snake segment \n including V-orbit centering and H-orbit compensation"    font "roman,13"   # offset 0,+.7    

 set xlabel "G{/Symbol g}"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set x2label "E [GeV]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "Spin rotation  [deg]"             font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" mirror      #offset 0,-.6

#set key above maxcol 2 1 
#set key tm c maxrows 1 width 4
set key t r maxrow 5   width -3
set key font "roman, 8"  samplen 1

set grid

pi = 4.*atan(1.)
deg = 180./pi
am = 938.27203
G = 1.79284735

Ggmi = 45 ; Ggma = 300e3/am*G

#set samples 10000
set xrange  [Ggmi:Ggma]
set x2range [Ggmi/G*am/1e3:Ggma/G*am/1e3]

#set xtics add (487 487)
#set arrow from 487,-180 to 487,-170 nohead

 plot \
  "zgoubi.SPNPRT.Out"  u (21==3 && 23==1552 ? $18 : 1/0):(acos($15)/3.1416*180.)  w lp pt 4 ps .4 lt 2 lc 1 lw 1 tit "SPINR"  ,\
  "zgoubi.SPNPRT.Out"  u (21==3 && 23==1552 ? $18/G*am/1e3 : 1/0):(acos($15)/3.1416*180.)  axes x2y1 w lp ps 0. lw 0. notit  

 set samples 10000
 set terminal postscript eps blacktext color enh size 8.3cm,5cm "Times-Roman" 12
 set output "gnuplot_SpinRotInRHICSegment.eps"
 replot 
 set terminal X11
 unset output

pause 2

 exit
