C So to get dp/p ranges for zgoubi simulations, 
C compute dp/p=dE/E/bta^2 for various E :  1.5 -> 25 GeV, given DgammaG. 
C Assuming DgG = .25 allows deriving value for DE = DgG x M/G

      am = 0.93827231e9
      G = 1.7928474
      DgG = 0.25

      DE = DgG * am / G
      write(*,*) ' '
      write(*,*) ' We consider full delta-Energie  DE of : ',DE/1e6,
     >' MeV,   whatever E.   dp/p follows : ' 
      write(*,*) ' '

      T = 1.5e9
      Tinc = .5e9
      do i = 1, 40
        E = T + am
        p = sqrt(E*E - am*am)
        bta = p / E   
        dpp = DE/E /(bta*bta)
        gG = E/am * g
        write(*,*) ' kin_E, dE/E, dp/p, bta, gamma*G : '
     >    ,T/1e6, DE/E, dpp ,bta,gG
        T = T + Tinc
      enddo

      stop
      end










