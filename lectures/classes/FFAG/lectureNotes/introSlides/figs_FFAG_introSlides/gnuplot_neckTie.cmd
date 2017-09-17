
set title "Necktie diagram, MARK I style of FFAG"    font "roman,16"   # offset 0,+.7    

 set xlabel "{/Symbol f}_1 [/2{/Symbol p}]"        font "roman,16"   # offset +4,-.5 rotate by +20  
 set ylabel "{/Symbol f}_2 [/2{/Symbol p}]"        font "roman,16"   # offset +4,-.5 rotate by +20  

 set xtics  font "roman,12" mirror
 set ytics  font "roman,12" mirror      #offset 0,-.6


plot [0:1] [0:1] 'fort.88' u ($1):($3==1 ? $2 : 1/0) w p ps .3 notit

 set terminal postscript eps blacktext color enh size 9.cm, 9.cm "Times-Roman" 12
 set output "gnuplot_neckTie.eps"
 replot
 set terminal X11
 unset output

      pause 2

 exit

