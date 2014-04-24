

#set xrange [.565:.585]

plot [0.:.5][] "spectrum.Out" w l 



pause 8
   
 set terminal postscript eps blacktext color
 set output "gnuplot_spectrum.eps"
 replot
 set terminal X11
 unset output

