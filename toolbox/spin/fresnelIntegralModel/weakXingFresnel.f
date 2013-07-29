      implicit double precision (a-h,o-z)
      parameter (pi = 4.d0 * atan(1.d0), pi2 = pi/2.d0)
      parameter(am = 938.27203d6, c= 2.99792458d8, G = 1.7928474)
      parameter (lunR=7, lunW=8)

      open(unit=lunR,file='weakXingFresnel.in')
      open(unit=lunW,file='weakXingFresnel.out')

      read(lunR,*) V        ! volts
      read(lunR,*) phis     ! rad
      read(lunR,*) wn2      ! Jn^2, resonance strength (width), can be from prior Froissard-Stora tracking
      read(lunR,*) n        ! such that gammaG=n+/-qz
      read(lunR,*) qz       ! algebraic value, such that gammaG=n+qz == n+/-abs(qz)
      read(lunR,*) miTurn
      read(lunR,*) maTurn

c      V = 6.d3        
c      phis = 0.4873884165731143248
c      wn2 =2.8635757d-9      
c      n  = 8           !  8-qz
c      qz = 3.62074437
c      miTurn = -1200
c      maTurn = 1400

      dW = V*sin(phis)
      a  = G * dW / am / (2.d0*pi)       ! =xing speed = Gdg/dtta
         write(lunW,*) '# Xing speed a = Gdg/dtta = ',a
         write(lunW,*) '# |J_n|^2 = ',wn2,',    pi|J_n|^2/a = ',pi*wn2/a

      gGr = dble(n) + qz
      Er = gGr /G *am
      Tr = Er-am
      dTdtta = a/G*am  *2.d0*pi  

      w2sa = pi * (wn2 /a) 
      rap = sqrt(a/pi) * (2.d0 *pi)
CCCCCCCCCCC VERIFIER fints !!*************
      do i = miTurn, 0, 1
        y = rap * dble(i)
        if(-y.gt.5.d0)   stop ' y too large' 
        s2p = w2sa * ((0.5d0-fintc(-y))**2 + (0.5d0-fints(-y))**2)
        T = Tr + dTdtta * dble(i)
        write(lunW,*) T/1d6 , sqrt(1-s2p), i
        write(88,*) T/1d6 , sqrt(1-s2p), i 
      enddo
      do i = 1, maTurn, 1
        y = rap * dble(i)
        if(y.gt.5.d0)   stop ' y too large' 
        s2p = w2sa * ((0.5d0+fintc(y))**2 + (0.5d0+fints(y))**2)
        T = Tr + dTdtta * dble(i)
        write(lunW,*) T/1d6 , sqrt(1-s2p), i
        write(88,*) T/1d6 , sqrt(1-s2p), i
      enddo               

      write(lunW,*) '# Resonant kin-E = ', Tr/1d6,' MeV'
      write(lunW,*) '# at theoretical turn# E_r/(dW/dTurn) = ', Tr/dW

      stop
      end
      function fintc(x)
      implicit double precision (a-h,o-z)
      parameter (pi = 4.d0 * atan(1.d0), pi2 = pi/2.d0)
      integer*8 nfctrl

      fintc = 0.d0
      f1 = fintc
      n = 0
 1    continue
        fac = x**(4*n+1) * (-1.d0)**n * pi2**(2*n) / dble(4*n+1) 
        m = 2
 2      continue
          if(m.gt.2*n) goto 3
          fac = fac / dble(m)
          m = m+1
          goto 2      
 3      continue
        fintc = fintc + fac
        if(abs(fintc-f1).le.1d-6*f1) goto 10
        f1 = fintc
        n = n+1
        goto 1

 10   return
      end
      function fints(x)
      implicit double precision (a-h,o-z)
      parameter (pi = 4.d0 * atan(1.d0), pi2 = pi/2.d0)
      integer*8 nfctrl

      fints = 0.d0
      f1 = fints
      n = 0
 1    continue
        fac = x**(4*n+3) * (-1.d0)**n * pi2**(2*n+1) / dble(4*n+3)
        m = 2
 2      continue
          if(m.gt.2*n+1) goto 3
          fac = fac / dble(m)
          m = m+1
          goto 2      
 3      continue
        fints = fints + fac
        if(abs(fints-f1).le.1d-6*f1) goto 10
        f1 = fints
        n = n+1
        goto 1

 10   return
      end
