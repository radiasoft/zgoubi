
      set title "POLARIZATION PROFILE"  font "roman,16" 

   set xlabel "\epsilon_y/sigma" font "roman,14" 
   set ylabel "<Sy>"             font "roman,14" offset 1,0 
   set y2label "#particles"      font "roman,14" 

      set xtics     mirror font "roman,11" 
      set ytics   nomirror font "roman,11" 
      set y2tics  nomirror font "roman,11" 

      set bmargin 4

      set key font "roman,15" bot left 
      set grid



# plot [][] \
#     "../231+Qy/profileS.out"  u ($1 < 40 ? $2 : 1/0):(abs($5))       w lp lw 1.5 ps .9  tit "231+Qy",\
#     "profileS.out"  u ($1 < 40 ? $2 : 1/0):(abs($5))       w lp  lw 1.5 ps .9  tit "411-Qy",\
#     "../393+Qy/profileS.out"  u ($1 < 40 ? $2 : 1/0):(abs($5))       w lp  lw 1.5 ps .9  tit "393+Qy"

 plot [][] \
     "profileS.out"  u ($2):(abs($5))                 w lp  lw 1.5 ps .9  tit "   Profile" ,\
     "profileS.out"  u ($2):($10)      axes x1y2      w lp  lw 1.5 ps .9  tit "   # prtcls" ,\
     "profileS.out"  u ($2):($12/10)      axes x1y2      w lp  lw 1.5 ps .9  tit "     Sum_prtcls/10"

       set terminal postscript eps blacktext color  enh size 8.3cm,4.4cm "Times-Roman" 14  
       set output "gnuplot_profileS.eps"  
       replot  
       set terminal X11  
       unset output  

pause 888

exit

   set y2label "#particles"      font "roman,14" 
      set y2tics  nomirror font "roman,11" 

#plot [][] \
#     "../231+Qy/profileS.out"  u ($1>1 && $1 < 40 ? $2 : 1/0):($9)  axes x2y2 w lp ps .5  notit   ,\
#     "profileS.out"  u ($1>1 && $1 < 40 ? $2 : 1/0):($10)  axes x2y2 w lp ps .5  notit   ,\
#     "../393+Qy/profileS.out"  u ($1>1 && $1 < 40 ? $2 : 1/0):($10)  axes x2y2 w lp ps .5  notit   


pause 888
