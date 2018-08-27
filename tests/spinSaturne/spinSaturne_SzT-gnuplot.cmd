set terminal postscript eps blacktext color enh "Times-Roman, 14"
set output "spinSaturne_SzT.eps"

set key maxcol 1
set key c l

#set logscale y

set tit "Spin component(s) vs. Turn number\n  (from zgoubi.fai) \n" font "Times-Roman, 14"

set xtics  nomirror font "Times-Roman, 14"
set x2tics nomirror font "Times-Roman, 14"
set ytics    mirror font "Times-Roman, 14"

set xlabel  'turn' font "Times-Roman, 16"
set x2label 'G{/Symbol g}' font "Times-Roman, 16"
set ylabel  'SZ' font "Times-Roman, 16"

# ELECTRON
# am = 0.5109989
# G = 1.159652e-3
# PROTON
am = 938.27208e6
G = 1.79284735

m2cm = 100.
MeV2eV = 1.e6
c = 2.99792458e8

V = 6000.
phs = 0.23635651
dEdt = V * sin(phs)
Ei = sqrt((5018.67e-3 * c) ** 2 + am ** 2)
Ef = Ei + 3000 * dEdt
Ggi = G * Ei / am
Ggf = G * Ef / am

set x2range[Ggi:Ggf]

plot for [i=1:3:1] \
   "zgoubi.fai" u ($26== i ? $38 : 1/0):($22) w p pt 7 ps .3 tit "SZ".i.""

#,\
#   "zgoubi.fai" u ($26== i ? G*(Ei + ($38-1)*dEdt)/am : 1/0):($22) axes x2y1 w p pt 2*i+1 ps .4 lc i  notit

exit

