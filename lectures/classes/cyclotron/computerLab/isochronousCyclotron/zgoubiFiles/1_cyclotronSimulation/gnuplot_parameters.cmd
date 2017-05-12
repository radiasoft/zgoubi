#plot "my.dat" every A:B:C:D:E:F
#A: line increment
#B: data block increment
#C: The first line        (line count starts at 0 !!)
#D: The first data block
#E: The last line
#F: The last data block

 set title "ISOCHRONOUS CYCLOTRON"    font "roman,14"   # offset 0,+.7    

 set xlabel "Radius [cm]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "TOF [{/Symbol m}s]"             font "roman,13"   #offset -0,-1 rotate by -20 
 set y2label "Kin. En. [MeV]"             font "roman,13"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" nomirror
 set x2tics  font "roman,12" nomirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics  font "roman,12" nomirror      #offset 0,-.6

set key t r maxrow 3   width -3
set key font "roman, 8"  samplen 1

set grid

pi = 4.*atan(1.)
deg = 180./pi
Ncell = 4

#set samples 10000
set xrange  []
set yrange []

 plot \
 'zgoubi.fai' u ($3):($15 *Ncell)       w l lc 1 tit 'TOF' ,\
 'zgoubi.fai' u ($3):($24) axes x1y2    w l lc 3 tit 'Kin. E' 

 set terminal postscript eps blacktext color enh size 8.3cm,5cm "Times-Roman" 12
 set output "gnuplot_parameters_TOF.eps"
 replot
 set terminal X11
 unset output

pause 2

 exit
