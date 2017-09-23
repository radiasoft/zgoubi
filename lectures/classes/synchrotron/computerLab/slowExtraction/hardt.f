      pi = 4.* atan(1.)

C At ESExtr
      D = 3.871
      Dp = -.617
      bet = 15.9
      alp = -0.351

      S =     13.95
      xi=-3.58  
      write(*,*) xi
C sextu to ESExtr
C      dmu = 228. /180. * pi
      dmu = 250. /180. * pi
C Separatrices at sextu, vol.I, p.19
      alpha = 180. /180. * pi 
 
      arg = alpha - dmu

      sqrb = sqrt(bet)

      Dn = D/sqrb
      Dpn = (alp*D + bet*Dp)/sqrb
     
      shl = Dn * cos(arg) + Dpn * sin(arg)
      shr = -4*pi*xi/S

      write(*,*) shl, shr

      stop 
      end

