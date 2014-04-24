
 set xtics font "roman,8"
 set ytics font "roman,14"

#set tit  'Fresnel integral model'  font "roman,25"
set xlab 'kin-E (MeV)'                    font "roman,25"
set ylab 'Sz'               font "roman,25"



plot  \
        'fort.88'  us ($1):($2)

pause 2

 set term post eps enh  size 8.3cm,4cm "helvetica,14" blacktext color
 set output "gnuplot_fresnelIntModel.eps"
 replot
 set terminal X11
 unset output


exit


