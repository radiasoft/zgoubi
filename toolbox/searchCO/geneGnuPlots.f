      character*120 cmmnd
      logical okGnu, makCO, makBta

C Write a gnuplot file for CO plot, and execute it
      okGnu= makCO()
      cmmnd = 'gnuplot < ./gnuplot_CO.cmd'
      write(*,*) '---------------------------------------'
      write(*,*) cmmnd
      call system(cmmnd)

      okGnu= makBta()
      cmmnd = 'gnuplot < ./gnuplot_betaFromMatrix.cmd'
      write(*,*) '---------------------------------------'
      write(*,*) cmmnd
      call system(cmmnd)

      stop
      end
      function makCO()
      logical makCO
      open(unit=8,file='gnuplot_CO.cmd')
      write(8,fmt='(a)') '# H and V CLOSED ORBITS '
      write(8,fmt='(a)') 'set xlabel "s (m)" font "roman,18"'
      write(8,fmt='(a)') 'set xtics font "roman,12"'
      write(8,fmt='(a)') 'set ytics font "roman,12"'
      write(8,fmt='(a)') ' '
      write(8,fmt='(a)') 'set grid'
      write(8,fmt='(a)') ' '
      write(8,fmt='(a)') 'set ylabel "x, y (m)" font "roman,18"'
      write(8,fmt='(2a)') 'set title "H & V closed orbits along ring"'
     > ,' font "roman,20"'
      write(8,fmt='(a)') ' '
C      write(8,fmt='(2a)') 'plot [:808][] \'  ! AGS
      write(8,fmt='(2a)') 'plot  \'
      write(8,fmt='(2a)') '"getDiffFromFai.out" '
     > //' u ($1 /100.):($2 /100.) w l lw 2  tit "x_{co}", \'
      write(8,fmt='(2a)') '"getDiffFromFai.out" '
     > //' u ($1 /100.):($3 /100.) w l lw 2  tit "y_{co}"'
      write(8,fmt='(a)') ' '
      write(8,fmt='(a)') 'pause 3'
      write(8,fmt='(a)') ' '
      write(8,fmt='(2a)') 'set terminal postscript eps blacktext color'
     > //' enh size 8.3cm,4cm "Times-Roman" 12 '
      write(8,fmt='(a)') ' set output "gnuplot_CO.eps" '
      write(8,fmt='(a)') ' replot '
      write(8,fmt='(a)') ' set terminal X11 '
      write(8,fmt='(a)') ' unset output '
      write(8,fmt='(a)') ' '

      write(8,fmt='(a)') ' '
      write(8,fmt='(a)') 'set ylabel "D_x, D_y (m)" font "roman,18"'
      write(8,fmt='(2a)') 'set title "H & V dispersion along ring"'
     > ,' font "roman,20"'
      write(8,fmt='(a)') ' '
C      write(8,fmt='(2a)') 'plot [:808][] \'  !AGS
      write(8,fmt='(2a)') 'plot  \'
      write(8,fmt='(2a)') '"getDiffFromFai.out" '
     > //' u ($1 /100.):($4 /100.) w l lw 2  tit "D_{x}",  \'
      write(8,fmt='(2a)') ' "getDiffFromFai.out" '
     > //' u ($1 /100.):($5 /100.) w l lw 2  tit "D_{y}"'
      write(8,fmt='(a)') ' '
      write(8,fmt='(a)') 'pause 3'
      write(8,fmt='(a)') ' '
      write(8,fmt='(2a)') 'set terminal postscript eps blacktext color'
     > //' enh size 8.3cm,4cm "Times-Roman" 12 '
      write(8,fmt='(a)') ' set output "gnuplot_DxDy.eps" '
      write(8,fmt='(a)') ' replot '
      write(8,fmt='(a)') ' set terminal X11 '
      write(8,fmt='(a)') ' unset output '
      write(8,fmt='(a)') ' '

      write(8,fmt='(a)') 'exit '
      close(8)
      makCO = .true.
      return
      end
      function makBta()
      logical makBta
      open(unit=8,file='gnuplot_betaFromMatrix.cmd')
  
      write(8,fmt='(a)') ' set xtics font "roman,18"'
      write(8,fmt='(a)') ' set ytics font "roman,18"'
      write(8,fmt='(a)') ' '
      write(8,fmt='(a)') 'set grid'
      write(8,fmt='(a)') ' '
      write(8,fmt='(a)') 'set tit  "beta_x,y versus s" font "roman,25"'
      write(8,fmt='(a)') 'set xlab "s (m)"             font "roman,25"'
      write(8,fmt='(2a)') 'set ylab "beta_x, -beta_y (m)" ',
     >' font "roman,25"'
C      write(8,fmt='(a)') 'plot [0:807.1][]\'   !AGS
      write(8,fmt='(a)') 'plot \'
      write(8,fmt='(2a)') ' "betaFromMatrix.out" ',
     >'us ($1/100.):($3) w l lw 2  tit "beta_x" , \'
      write(8,fmt='(2a)') ' "betaFromMatrix.out" ',
     >'us ($1/100.):($7 * (-1)) w l lw 2  tit "-beta_y" '
      write(8,fmt='(a)') ' '
      write(8,fmt='(2a)') ' set term post eps enh  size 8.3cm,4cm ',
     >'"helvetica,14" blacktext color'
      write(8,fmt='(a)') ' set terminal postscript eps blacktext color'
      write(8,fmt='(a)') ' set output "gnuplot_betxy.eps"'
      write(8,fmt='(a)') ' replot'
      write(8,fmt='(a)') ' set terminal X11'
      write(8,fmt='(a)') ' unset output'
      write(8,fmt='(a)') ' '
      write(8,fmt='(a)') ' '
      write(8,fmt='(a)') 'set tit  "D_x, D_y versus s" font "roman,25"'
      write(8,fmt='(a)') 'set xlab "s (m)"             font "roman,25"'
      write(8,fmt='(a)') 'set ylab "D_x, D_y (m)"      font "roman,25"'
C      write(8,fmt='(a)') 'plot [0:807.1][-2.5:1.]\'  !AGS
      write(8,fmt='(a)') 'plot \'
      write(8,fmt='(2a)') ' "betaFromMatrix.out" us ($1/100.):($4)',
     >' w l lw 2  tit "D_x" , \'
      write(8,fmt='(2a)') ' "betaFromMatrix.out" us ($1/100.):($8)',
     >' w l lw 2  tit "D_y" '
      write(8,fmt='(a)') ' '
      write(8,fmt='(2a)') ' set term post eps enh  size 8.3cm,4cm ',
     >'"helvetica,14" blacktext color'
      write(8,fmt='(a)') ' set terminal postscript eps blacktext color'
      write(8,fmt='(a)') ' set output "gnuplot_Dxy.eps"'
      write(8,fmt='(a)') ' replot'
      write(8,fmt='(a)') ' set terminal X11'
      write(8,fmt='(a)') ' unset output'
      write(8,fmt='(a)') ' '
      write(8,fmt='(a)') ' '
      write(8,fmt='(a)') ' '
      write(8,fmt='(2a)') 'set tit  "Phase advance [2pi] versus s " ',
     >' font "roman,25"'
      write(8,fmt='(a)') 'set xlab "s (m)"         font "roman,25"'
      write(8,fmt='(2a)') 'set ylab "mu_x/2pi, mu_y/2pi (m)" ',
     >'font "roman,25"'
C      write(8,fmt='(a)') 'plot [0: 12*  807.1/12][]\'   !AGS
      write(8,fmt='(a)') 'plot \'
      write(8,fmt='(2a)') ' "betaFromMatrix.out"       ',
     >'us ($1/100.):($10) w l lw 2  tit "mu_x" , \'
      write(8,fmt='(2a)') ' "betaFromMatrix.out"       ',
     >'us ($1/100.):($11) w l lw 2  tit "mu_y" '
      write(8,fmt='(a)') ' '
      write(8,fmt='(2a)') ' set term post eps enh  size 8.3cm,4cm',
     >' "helvetica,14" blacktext color'
      write(8,fmt='(a)') ' set terminal postscript eps blacktext color'
      write(8,fmt='(a)') ' set output "gnuplot_muxy.eps"'
      write(8,fmt='(a)') ' replot'
      write(8,fmt='(a)') ' set terminal X11'
      write(8,fmt='(a)') ' unset output'
      write(8,fmt='(a)') ' '
      write(8,fmt='(a)') 'exit'
      write(8,fmt='(a)') ' '

      close(8)
      makBta = .true.
      return
      end
