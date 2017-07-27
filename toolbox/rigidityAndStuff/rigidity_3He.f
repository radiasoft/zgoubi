      implicit double precision (a-h,o-z)
      parameter (c=2.99792458d8) 

      am = 2808.39148d6
      Q = 2.
      N = 3  ! nucleons
      G = -4.18415382
      write(6,*) '  Mass, Q, # Nucs : ',am, Q, N

 16   continue
      write(6,*) '  Give Ggamma : '
      read(*,*,err=16) Gg
           
         Gg = - abs(Gg)

      gamma = Gg/G
      E = gamma * am
      p = sqrt(E*E - am*am)
      brho = p / Q / c
      Ek = E - am
      bta = p / sqrt(p*p + am*am) 
            
      write(6,*) 'Rigidity, gamma : ',brho,gamma
      write(6,*) 'Energie, total, kinetic (GeV) : ',E/1d9,Ek/1d9
      write(6,*) 'Energie/N, total, kinetic (GeV/N) : ',E/1d9/N,Ek/1d9/N
      write(6,*) 'beta, beta*gamma, Ggamma : ',bta, bta*E/am, G*E/am
      write(6,*) ' '
      write(6,*) ' '
      write(6,*) ' '

        stop

      am = 2808.39148d6
      Q = 2.
      N = 3  ! nucleons
      G = -4.1841538
      write(6,*) '  Mass, Q, # Nucs : ',am, Q, N

 14   continue
      write(6,*) '  Give rigidity : '
      read(*,*,err=14) brho

      p = Q * brho * c
      E = sqrt(p*p + am*am)
      Ek = E - am
            
      write(6,*) 'Rigidity : ',brho
      write(6,*) 'Energie, total, kinetic (GeV) : ',E/1d9,Ek/1d9
      write(6,*) 'Energie/N, total, kinetic (GeV/N) : ',E/1d9/N,Ek/1d9/N
      write(6,*) 'beta, beta*gamma, Ggamma : ',bta, bta*E/am, G*E/am
      write(6,*) ' '
      write(6,*) ' '
      write(6,*) ' '

 15   continue
      write(6,*) '  Give Ggamma : '
      read(*,*,err=15) Gg

      gamma = Gg/G
      E = gamma * am
      p = sqrt(E*E - am*am)
      brho = p / Q / c
      Ek = E - am
            
      write(6,*) 'Rigidity : ',brho
      write(6,*) 'Energie, total, kinetic (GeV) : ',E/1d9,Ek/1d9
      write(6,*) 'Energie/N, total, kinetic (GeV/N) : ',E/1d9/N,Ek/1d9/N
      write(6,*) 'beta, beta*gamma, Ggamma : ',bta, bta*E/am, G*E/am
      write(6,*) ' '
      write(6,*) ' '
      write(6,*) ' '

      stop
      end

 
