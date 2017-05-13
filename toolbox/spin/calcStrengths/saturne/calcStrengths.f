      parameter (pi = 4.*atan(1.), am = 938.27203e6, deg2rd = pi/180.)
      parameter (G = 1.7928474)

      circ = 105.557524    ! m
      Vp = 6000.           ! V
      phis = 27.92529925320347564 * deg2rd  ! rad
      Bdip =  0.157776     ! T
      Brho = 1.

      DE = Vp * sin(phis)      
      alpha = G * DE / am / (2. * pi)
      rho = 1./Bdip
      R = circ / (2.*pi)
      Bdot = DE/(2.*pi*R*rho)

      write(*,*) ' Present conditions : '
      write(*,*) '    xing speed alpha = dgamma/dtta =   ',alpha
      write(*,*) '    dB/dt                         =    ',Bdot

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 8-nu_z '
      write(*,*) ' ' 
      ez = 0.01e-6
      pini = 1.
      pfin = 0.9895
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ez,pini,pfin,A2,A2/ez,A2*2.*alpha/pi,A2*2.*alpha/pi/ez


      stop
      end





