
 set title "Polarization, eRHIC storage ring \n From polSpectrumFromFai.out"    font "roman,16"   # offset 0,+.7    

 set xlabel "a{/Symbol g}"        font "roman,18"   # offset +4,-.5 rotate by +20  
 set ylabel "S_l,  turn survival (rel.)"          font "roman,18"   #offset -0,-1 rotate by -20 
 set y2label "{/Symbol s}_{/Symbol d}E/E_{ref}" font "roman,18"   #offset -0,-1 rotate by -20 

 set xtics  font "roman,12" mirror
 set ytics  font "roman,12" nomirror      #offset 0,-.6
 set y2tics  font "roman,12" nomirror      #offset 0,-.6

set key t l
set key maxrow 1 # width 1
set key font "roman, 10"  # samplen 1
#set key  spacin 2

set grid
set samples 10000 

betx = 0.63309
bety = 0.07606
cm2m=0.01
mrd2rd=0.001
m2mum = 1e6
gma = 18e3/.511


# set xrange [21.8:23.2] 
 set xrange [39.8:41.2] 
set yrange [-1:1.3]

aGammaRef = 40.5
sigOv2 = .5    # In order for the xerror bar width to be 1 sigma of a.gamma (not +\- 1sigma)
nbPass = 225000

 plot \
     "polSpectrumFromFai.out" u ($5):($6) with l linecolor rgb "blue" tit 'S_l' ,\
     "polSpectrumFromFai.out" u ($5):($10/nbPass) with l linecolor rgb "red" tit 'turn srvvl' ,\
     "polSpectrumFromFai.out" u ($5):($3) axes x1y2 w l linecolor rgb "green" tit '{/Symbol s}_{/Symbol d}E/E_{ref}'

# plot \
#     "polSpectrumFromFai.out" u ($2):($6):($3 *sigOv2) with xerrorbars linecolor rgb "blue" notit ,\
#     "polSpectrumFromFai.out" u ($2):($6) with l linecolor rgb "blue" notit ,\
#     "polSpectrumFromFai.out" u ($2):($3/$5) axes x1y2 w l linecolor rgb "red" tit "{/Symbol s}_{{/Symbol d}E}/E_{ref}" 

 set terminal postscript eps blacktext color enh "Times-Roman" 12
 set output "gnuplot_polSpec.eps"
 replot
 set terminal X11
 unset output
  
pause 20

  

