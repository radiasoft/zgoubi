      parameter (pi = 4.*atan(1.), am = 938.27203e6, deg2rd = pi/180.)
      parameter (G = 1.7928474)

      circ =  3833.845900000    ! m
      Vp = 300.e3               ! kV
      phis = 20. * deg2rd       ! rad
      Bdip =  4.1230212E-03     ! T   for Brho=1T.m 
      Brho = 1.

      DE = Vp * sin(phis)      
      alpha = G * DE / am / (2. * pi)
      rho = Brho/Bdip
      R = circ / (2.*pi)
      Bdot = DE/(2.*pi*R*rho)

      write(*,*) ' Present conditions : '
      write(*,*) '    xing speed alpha = dgamma/dtta =   ',alpha
      write(*,*) '    dB/dt                         =    ',Bdot

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 318 - nu_z  (nu_z~0.6704))'
      write(*,*) ' ' 
      ez =  4.4e-8
      pini = 0.9838
      pfin = 0.9054
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin,
     >         '  p_final/p_final = ',pfin/pini
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.


      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 318 + nu_z  (nu_z~0.6704))'
      write(*,*) ' ' 
      ez =  6.28e-11
      pini = 0.9996
      pfin = 0.5461
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin,
     >         '  p_final/p_final = ',pfin/pini
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.



      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 330 - nu_z  (nu_z~0.6704))'
      write(*,*) ' ' 
      ez =  2.90e-11
      pini = 0.9998
      pfin = 0.5499
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin,
     >         '  p_final/p_final = ',pfin/pini
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.


      stop
      end





