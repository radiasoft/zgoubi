reset

 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set terminal postscript eps blacktext color
 set output "gnuplot_betaFctns.eps"
 
#------- Siziing margins and plots
NX=1; NY=3
###DX=0.15; DY=0.08; SX=0.7; SY=0.26
DX=0.1; DY=0.08; SX=0.7; SY=0.26
set bmargin 20*DY; set tmargin 20*DY; set lmargin 20.*DX; set rmargin 80*DX
## set the margin of each side of a plot as small as possible
## to connect each plot without space
set size SX*NX+DX*1.5,SY*NY+DY*1.8
#--------------------------------
 
set grid
set tit  'RHIC functions from Zgoubi'
xmin = 0.
xmax = 3860.
set xrange [xmin:xmax]

unset xtics
unset xlabel

#----------------------------------
set multiplot layout 3,1  rowsfirst
#----------------------------------

# TOP PLOT----------------------

set size SX,SY
set origin DX,DY+2*SY 

set ytics font "roman,10"  nomirror 
set y2tics font "roman,10"  nomirror 

set key default
set key top right spacing .75 box 2  
set label 1 '(a)' at graph 0.02,0.9 font ",14"
set ylab  "\beta_x (m)"  font "roman,14" tc lt 1 
set y2lab "\beta_y (m)"  font "roman,14" tc lt 3

cm2m = 0.01

plot [][] \
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin ? $2 : 1/0)           w l lc 1 tit "\beta_x"  ,\
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin ? $2 : 1/0)           w l lc 1 notit          ,\
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin && stringcolumn(1) eq 'DPUE' ? $2 : 1/0)           w p pt 6 ps .4 lc 1 notit ,\
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin ? $4 : 1/0) axes x1y2 w l lc 3 tit "\beta_y"  ,\
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin ? $4 : 1/0) axes x1y2 w l lc 3 notit  ,\
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin && stringcolumn(1) eq 'DPUE' ? $4 : 1/0) axes x1y2 w p pt 6 ps .4 lc 3 notit 


unset title
 

# MIDDLE PLOT----------------------
set size SX,SY
set origin DX,DY+SY

set ytics  font "roman,10"  nomirror 
set y2tics font "roman,10"  nomirror 

set key top right box 2 
set label 1 '(b)' at graph 0.02,0.9 font ',14'
set ylab  "Dx (m)"  font "roman,14" 
set y2lab "Dy (m)"  font "roman,14" 

plot \
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin ? $7 : 1/0)             w l lc 1  tit "Dx"  ,\
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin ? $7 : 1/0) axes x1y2   w l lc 1  notit     ,\
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin && stringcolumn(1) eq 'DPUE' ? $7 : 1/0)           w p pt 6 ps .4 lc 1 notit ,\
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin ? $9 : 1/0)             w l lc 3  tit "Dy"  ,\
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin ? $9 : 1/0) axes x1y2   w l lc 3  notit ,\
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin && stringcolumn(1) eq 'DPUE' ? $9 : 1/0) axes x1y2 w p pt 6 ps .4 lc 3 notit 
 

# BOTTOM PLOT----------------------
set size SX,SY
set origin DX,DY 

unset key
set key bot right box 2 
set xlab "s (m)"  font "roman,16" offset -5, 0.5
set xtics font "roman,10"  mirror 

set label 1 '(c)' at graph 0.02,0.9 font ',14'
set label '* Qx= 0.8257 Qy= 0.9445 Chrox= -24.13 Chroy= -90.11' at screen 0.1,0.025 font ',14'

set ylab  "x (m)"  font "roman,14"
set y2lab "y (m)"  font "roman,14"

plot \
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin ? $15 : 1/0)           w l lc 1  tit "x"  ,\
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin && stringcolumn(1) eq 'DPUE' ? $15 : 1/0)           w p pt 6 ps .4 lc 1 notit ,\
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin ? $17 : 1/0) axes x1y2 w l lc 3  tit "y"  ,\
 "zgoubi.TWISS.out" us ($13 *cm2m):($13 *cm2m<xmax && $13 *cm2m>xmin && stringcolumn(1) eq 'DPUE' ? $17 : 1/0) axes x1y2 w p pt 6 ps .4 lc 3 notit 

unset multiplot

system'okular gnuplot_opticalFctns.eps'

exit



 
