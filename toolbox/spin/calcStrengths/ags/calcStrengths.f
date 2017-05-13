      parameter (pi = 4.*atan(1.), am = 938.27203e6, deg2rd = pi/180.)
      parameter (G = 1.7928474)

      circ = 807.04378     ! m
      Vp = 290.e3          ! kV
      phis = 30. * deg2rd  ! rad
      Bdip =  0.01171255    
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
      write(*,*) ' gamma G = nu_z '
      write(*,*) ' ' 
      ez = 0.002e-6
      pini = 1.
      pfin = 0.98582
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
     >              ,' $ pi.|J_n|^2/alpha =',A2*2.,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = nu_z '
      write(*,*) ' ' 
      ez = 0.01e-6
      pini = 1.
      pfin = 0.9298
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
     >              ,' $ pi.|J_n|^2/alpha =',A2*2.,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = nu_z '
      write(*,*) ' ' 
      ez = 0.05e-6
      pini = 1. 
      pfin = 0.6725
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = nu_z '
      write(*,*) ' ' 
      ez = 0.2e-6
      pini = 0.9999 
      pfin = -5.22e-2
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 24 - nu_z '
      write(*,*) ' ' 
      ez =  0.1e-6
      pini = 1. 
      pfin = 0.997900 
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.



      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 12 + nu_z '
      write(*,*) ' ' 
      ez =  0.002e-6
      pini = 1. 
      pfin = 0.994876 
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.


      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 12 + nu_z '
      write(*,*) ' ' 
      ez =  0.1e-6
      pini = 1. 
      pfin = 0.75800 
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.


      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 12 + nu_z '
      write(*,*) ' ' 
      ez =  2e-6
      pini = 1. 
      pfin = -0.842770 
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.


      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 24 + nu_z '
      write(*,*) ' ' 
      ez =  2.e-6
      pini = 0.999
      pfin = 0.857
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 24 + nu_z '
      write(*,*) ' ' 
      ez =  30.e-6
      pini = 0.987
      pfin = -0.269
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 36 - nu_z '
      write(*,*) ' ' 
      ez =  0.05e-6
      pini = 0.99983
      pfin = .1545
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 36 + nu_z '
      write(*,*) ' ' 
      ez = 0.0001e-6
      pini = 1. 
      pfin = 0.98682
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
     >              ,' $ pi.|J_n|^2/alpha =',A2*2.,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 36 + nu_z '
      write(*,*) ' ' 
      ez = 0.002e-6
      pini = 1. 
      pfin = 0.751
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 36 + nu_z '
      write(*,*) ' ' 
      ez = 0.02e-6
      pini = 1
      pfin = -0.522
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 48 - nu_z, e_z/pi=.125 '
      write(*,*) ' ' 
      ez = .125e-6
      pini = 1.
      pfin =  0.9095 
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 48 - nu_z, e_z/pi=.25 '
      write(*,*) ' ' 
      ez =  .25e-6
      pini =  0.9998
      pfin =    0.812 
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 48 - nu_z, e_z/pi=.5 '
      write(*,*) ' ' 
      ez =  .5e-6
      pini =   0.9996
      pfin =   0.6515  
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 48 - nu_z, e_z/pi=1.e-6 '
      write(*,*) ' ' 
      ez =  1e-6
      pini =  0.9993 
      pfin =   0.3545 
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 48 - nu_z, e_z/pi=2e-6 '
      write(*,*) ' ' 
      ez =  2e-6
      pini = 0.9985 
      pfin =  -0.1078  
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' e_z/pi = ', ez
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/e_z   =',A2/ez,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/ez =',A2*2.*alpha/pi/ez,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Integer resonances

      write(*,*) 
      write(*,*) '----------------------- '
      write(*,*) ' gamma G = 9,  z = 0.093mm '
      write(*,*) ' ' 
      z = 0.093e-3
      z2 = z*z
      pini = 1.
      pfin = .997457
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' z^2 = ', z2
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/z^2   =',A2/z2,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/z^2 =',A2*2.*alpha/pi/z2,' $'
     >              ,' $ pi.|J_n|^2/alpha =',A2*2.,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '----------------------- '
      write(*,*) ' gamma G = 9,  z = 0.93mm '
      write(*,*) ' ' 
      z = 0.93e-3
      z2 = z*z
      pini = 1.
      pfin = .762110
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' z^2 = ', z2
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/z^2   =',A2/z2,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/z^2 =',A2*2.*alpha/pi/z2,' $'
     >              ,' $ pi.|J_n|^2/alpha =',A2*2.,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '----------------------- '
      write(*,*) ' gamma G = 9,  z = 9.3mm '
      write(*,*) ' ' 
      z = 9.3e-3
      z2 = z*z
      pini = .992876
      pfin = -.995
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' z^2 = ', z2
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/z^2   =',A2/z2,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/z^2 =',A2*2.*alpha/pi/z2,' $'
     >              ,' $ pi.|J_n|^2/alpha =',A2*2.,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '----------------------- '
      write(*,*) ' gamma G = 9,  z = 1.385 '
      write(*,*) ' ' 
      z = 1.385e-3
      z2 = z*z
      pini = 1.
      pfin = 0.9941
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' z^2 = ', z2
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/z^2   =',A2/z2,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/z^2 =',A2*2.*alpha/pi/z2,' $'
     >              ,' $ pi.|J_n|^2/alpha =',A2*2.,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '----------------------- '
      write(*,*) ' gamma G = 9,  z = 2.77 '
      write(*,*) ' ' 
      z = 2.77e-3
      z2 = z*z
      pini = 0.9992
      pfin = -0.412
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' z^2 = ', z2
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/z^2   =',A2/z2,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/z^2 =',A2*2.*alpha/pi/z2,' $'
     >              ,' $ pi.|J_n|^2/alpha =',A2*2.,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 12 '
      write(*,*) ' ' 
      z = 13.84e-3
      z2 = z*z
      pini = 0.9999
      pfin = 5.674e-2
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' z^2 = ', z2
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/z^2   =',A2/z2,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/z^2 =',A2*2.*alpha/pi/z2,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 13 '
      write(*,*) ' ' 
      z = 13.84e-3
      z2 = z*z
      pini = 0.9999
      pfin = 0.15745
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' z^2 = ', z2
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/z^2   =',A2/z2,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/z^2 =',A2*2.*alpha/pi/z2,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 23 '
      write(*,*) ' ' 
      z = 1.36e-2
      z2 = z*z
      pini = 0.99975
      pfin = -0.6468
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' z^2 = ', z2
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/z^2   =',A2/z2,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/z^2 =',A2*2.*alpha/pi/z2,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 27 '
      write(*,*) ' ' 
      z = 1.384e-3
      z2 = z*z
      pini = 0.9998
      pfin = -0.167355 
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' z^2 = ', z2
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/z^2   =',A2/z2,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/z^2 =',A2*2.*alpha/pi/z2,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 45 '
      write(*,*) ' ' 
      z = 0.0275e-3
      z2 = z*z
      pini = 1.
      pfin = 0.99596
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' z^2 = ', z2
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/z^2   =',A2/z2,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/z^2 =',A2*2.*alpha/pi/z2,' $'
     >              ,' $ pi.|J_n|^2/alpha =',A2*2.,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      write(*,*) 
      write(*,*) '-------------------'
      write(*,*) ' gamma G = 45 '
      write(*,*) ' ' 
      z = 0.275e-3
      z2 = z*z
      pini = 0.999755
      pfin = 0.62962
      pfi1 = pfin/pini + 1.
      A2 = -alog(pfi1/2.)
      write(*,*) ' z^2 = ', z2
      write(*,*) ' p_init = ',pini,'   p_final = ',pfin
      write(*,*) '$ A^2   =',A2,' $','$ A^2/z^2   =',A2/z2,' $'
      write(*,*) ' $ |J_n|^2 =',A2*2.*alpha/pi,' $'
     >              ,' $ |J_n|^2/z^2 =',A2*2.*alpha/pi/z2,' $'
     >              ,' $ pi.|J_n|^2/alpha =',A2*2.,' $'
      write(*,*) ' hence pf/pi = ',2.*exp(-A2) - 1.

      stop
      end





