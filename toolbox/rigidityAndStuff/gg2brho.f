      implicit double precision (a-h,o-z)

      am = 938.27203d6
      c = 2.99792458d8
      G = 1.7928474d0

      Gga = 4.5
      ga = Gga /  G
      T = am * (ga-1.d0)
      brho = sqrt(T*(T+2.d0*am))/c
      write(*,*) ' G.gamma, Brho : ',Gga,brho

      Gga = 45.5
      ga = Gga /  G
      T = am * (ga-1.d0)
      brho = sqrt(T*(T+2.d0*am))/c
      write(*,*) ' G.gamma, Brho : ',Gga,brho

      Gga = 43.5
      ga = Gga /  G
      T = am * (ga-1.d0)
      brho = sqrt(T*(T+2.d0*am))/c
      write(*,*) ' G.gamma, Brho : ',Gga,brho

      brho = 75.8725752916200  
      p = brho *c
      E = sqrt(p*p + am*am)
      Gga = G * E/am       
      write(*,*) ' Brho, G.gamma : ',Gga,brho

      Gga = 377.57
      ga = Gga /  G
      T = am * (ga-1.d0)
      brho = sqrt(T*(T+2.d0*am))/c
      write(*,*) ' G.gamma, Brho : ',Gga,brho

      Qy = 29.6875
      GgaR = 411.d0 - Qy
      GgaIni = GgaR - dble(6)
      gaIni = GgaIni /  G
      T = am * (gaIni-1.d0)
      brhoIni = sqrt(T*(T+2.d0*am))/c
      write(*,*) ' G.gamma, Brho : ',GgaIni,brhoIni

      stop
      end
 

