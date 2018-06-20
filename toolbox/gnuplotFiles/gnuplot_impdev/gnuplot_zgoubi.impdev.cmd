
Brho = 3.3658103258583867e3 *1e-3  # T.m

      set xtics     mirror font "roman,11" 
      set ytics     mirror font "roman,11" 

      set key font "roman,15" t l
      set key maxcol 1
      set grid


      set title "Field and derivatives across magnet  \n BZ  "  font "roman,16"  


      unset ytics 
      set ytics    mirror font "roman,11" 
      set  nogrid

set xlabel "distance (cm)" font "roman,14" 
set ylabel "B_y (T)" font "roman,14" 

plot [-5:55] \
     'zgoubi.impdev.out' u 7:($10/10) w lp ps .04 notit

       set terminal postscript eps blacktext color  enh size 8.3cm,4.4cm "Times-Roman" 14  
       set output "gnuplot_impdev_BZ.eps"  
       replot  
       set terminal X11  
       unset output  

      unset ytics 
      set ytics     mirror font "roman,11" 
      unset y2tics
      unset y2label 

pause .5

      set title "Field and derivatives across magnet"  font "roman,16" 

set ylabel "dB_X/dX" font "roman,14" 

plot 'zgoubi.impdev.out' u 7:11 w lp ps .4 tit 'dBX/dX'

pause .5

      set title "Field and derivatives across magnet \n dBY/dX  "  font "roman,16" 

plot \
   'zgoubi.impdev.out' u ($410 ==1 ? $7 : 2/0):12 w lp ps .03 lc 1 tit 'dBY/dX 40.5MeV' ,\
   'zgoubi.impdev.out' u ($410 ==2 ? $7 : 2/0):12 w lp ps .03 lc 3 tit '       42' 

       set terminal postscript eps blacktext color  enh size 8.3cm,4.4cm "Times-Roman" 14  
       set output "gnuplot_impdev_dBYdX.eps"  
       replot  
       set terminal X11  
       unset output  

pause .5

      set title "Field and derivatives across magnet \n dBZ/dX  "  font "roman,16" 

plot \
   'zgoubi.impdev.out' u ($410 ==1 ? $7 : 2/0):13 w lp ps .03 lc 1 tit 'dBZ/dX 40.5MeV' ,\
   'zgoubi.impdev.out' u ($410 ==2 ? $7 : 2/0):13 w lp ps .03 lc 3 tit '       42' 

       set terminal postscript eps blacktext color  enh size 8.3cm,4.4cm "Times-Roman" 14  
       set output "gnuplot_impdev_dBZdX.eps"  
       replot  
       set terminal X11  
       unset output  

pause .5

      set title "Field and derivatives across magnet \n dBX/dY  "  font "roman,16" 

plot \
   'zgoubi.impdev.out' u ($410 ==1 ? $7 : 2/0):14 w lp ps .03 lc 1 tit 'dBX/dY 40.5MeV' ,\
   'zgoubi.impdev.out' u ($410 ==2 ? $7 : 2/0):14 w lp ps .03 lc 3 tit '       42' 

       set terminal postscript eps blacktext color  enh size 8.3cm,4.4cm "Times-Roman" 14  
       set output "gnuplot_impdev_dBXdY.eps"  
       replot  
       set terminal X11  
       unset output  


pause .5

      set title "Field and derivatives across magnet   "  font "roman,16" 

set xlabel "distance (cm)" font "roman,14" 
set ylabel "dB_Y/dY " font "roman,14" 
plot 'zgoubi.impdev.out' u 7:15 w lp ps .4 tit 'dBY/dY'

pause .5


set xlabel "distance (m)" font "roman,14" 
set ylabel "dBZ/dY /B{/Symbol r} [m^{-2}]" font "roman,14" 
plot 'zgoubi.impdev.out' u ($7/100.):($16 /Brho) w lp ps .4 tit 'dBZ/dY'

set terminal postscript eps blacktext color  enh  "Times-Roman" 14  
       set output "gnuplot_impdev_dBZdY.eps"  
       replot  
       set terminal X11  
       unset output  

pause .5

plot 'zgoubi.impdev.out' u 7:17 w lp ps .4 tit 'dBX/dZ'

pause .5

set xlabel "distance (m)" font "roman,14" 
set ylabel "dBY/dZ /B{/Symbol r} [m^{-2}]" font "roman,14" 
plot 'zgoubi.impdev.out' u ($7/100.):($18 /Brho) w lp ps .4 tit 'dBY/dZ'

set terminal postscript eps blacktext color  enh  "Times-Roman" 14  
       set output "gnuplot_impdev_dBYdZ.eps"  
       replot  
       set terminal X11  
       unset output  


pause .5

plot 'zgoubi.impdev.out' u 7:19 w lp ps .4 tit 'dBZ/dZ'

pause .5


plot 'zgoubi.impdev.out' u 7:20 w lp ps .4 tit 'd2BX/dX2'

pause .5

plot 'zgoubi.impdev.out' u 7:21 w lp ps .4 tit 'd2BY/dX2'

pause .5

      set title "Field and derivatives across magnet \n d2BZ/dX2  "  font "roman,16" 

plot \
   'zgoubi.impdev.out' u ($410 ==1 ? $7 : 2/0):22 w lp ps .03 lc 1 tit 'd2BZ/dX2 40.5MeV' ,\
   'zgoubi.impdev.out' u ($410 ==2 ? $7 : 2/0):22 w lp ps .03 lc 3 tit '         42'

       set terminal postscript eps blacktext color  enh size 8.3cm,4.4cm "Times-Roman" 14  
       set output "gnuplot_impdev_d2BZdX2.eps"  
       replot  
       set terminal X11  
       unset output  

pause .5




exit
