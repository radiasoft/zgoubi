      implicit double precision (a-h,o-z)
      parameter (am=938.27203D6, c=2.99792458D8, q=1.602176487D-19)
      parameter (G = 1.7928474)

! AGS 
      Et = am * 4.5/G
      T = Et - am
      br = sqrt(T * (T + 2.d0 * am))/c
      write(*,*) ' AGS, gG = 4.5 => T, Brho = ',T,br
      write(*,*) 

      Et = am * 5/G
      T = Et - am
      br = sqrt(T * (T + 2.d0 * am))/c
      write(*,*) ' AGS, gG = 5 => T, Brho = ',T,br
      write(*,*) 

      TInj = 1.5d9
      brInj = sqrt(TInj * (TInj + 2.d0 * am))/c
      write(*,*) ' AGS , T_Inj, Brho_Inj : ',TInj, brInj
      write(*,*) 

      TXtr = 24.3d9
      brXtr = sqrt(TXtr * (TXtr + 2.d0 * am))/c
      write(*,*) ' AGS , T_extr, Brho_extr : ',TXtr, brXtr
      write(*,*) ' AGS , Brho_extr/BrhoInj : ',brXtr/brInj
      write(*,*) 

      T = 23.8124d9
      br = sqrt(T * (T + 2.d0 * am))/c
      write(*,*) ' AGS , T, Brho : ',T, br
      write(*,*) 

      write(*,*) 
      write(*,*) 
      Et = am * 43.5/G
      T = Et - am
      br = sqrt(T * (T + 2.d0 * am))/c
      write(*,*) ' AGS, gG = 43.5 => T, Brho = ',T,br
      write(*,*) 

      Ttemp = T
      write(*,*) 
      Et = am * 46.5/G
      T = Et - am
      br = sqrt(T * (T + 2.d0 * am))/c
      write(*,*) ' AGS, gG = 46.5 => T, Brho = ',T,br
      write(*,*) 
      dEdn = 290000.d0 * sin(2.617993877991494365)
      write(*,*) ' dE/dn, number of turns 43.5/G -> 46.5 : ',
     >                                             dEdn,(T-Ttemp)/dEdn
      write(*,*) 
      write(*,*) 

! RHIC     
      Et = 144.d0 *am
      dE = (151.d0 -144.d0)*am
      xturn = dE / 150d3
      T = Et - am
      br = sqrt(T * (T + 2.d0 * am))/c
      write(*,*) ' RHIC, gamma = 144.  am, Brho, #turns : ',am,br,xturn
      write(*,*) 

      Et = 45.5/G*am
      T = Et - am
      br = sqrt(T * (T + 2.d0 * am))/c
      write(*,*) ' RHIC, injection, gG = 45.5 => T, Brho = ',T,br
      write(*,*) 

      Ttemp = T

      Et = 45.68/G*am
      T = Et - am
      br = sqrt(T * (T + 2.d0 * am))/c
      write(*,*) ' RHIC, gG=16+Qy => T, Brho = ',T,br
      write(*,*) 

      dEdn = 300000.d0 * sin(2.7925)
      write(*,*) ' dE/dn, number of turns 45.5/G -> 16+Qy : ',
     >                                             dEdn,(T-Ttemp)/dEdn
      write(*,*) 
      Et = 422/G*am
      T = Et - am
      br = sqrt(T * (T + 2.d0 * am))/c
      write(*,*) ' RHIC, gG=16+Qy => T, Brho = ',T,br
      write(*,*) 

      dEdn = 300000.d0 * sin(2.7925)
      write(*,*) ' dE/dn, number of turns 45.5/G -> 422/G : ',
     >                                             dEdn,(T-Ttemp)/dEdn
      write(*,*) 

      Et = am * 422./G
      T = Et - am
      br = sqrt(T * (T + 2.d0 * am))/c
      write(*,*) ' RHIC, top E, gG = 422 => T, Brho = ',T,br
      write(*,*) 


      stop
      end

