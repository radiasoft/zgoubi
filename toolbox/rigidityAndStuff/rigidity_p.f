      implicit double precision (a-h,o-z)
      parameter (c=2.99792458d8) 

      am = 938.27203d6 
      nq = 1
      G = 1.79284735
      write(6,*) ' p mass, #q : ', am , nq

 2    continue
      write(6,*) '  Give kinetic E (eV) : '
      read(*,*,err=2) Ek

      p = sqrt(Ek*(Ek + 2.d0*am))
      brho = p/dble(nq)/c
      bta = p / sqrt(p*p + am*am)
      E = Ek + am
            
      write(6,*) ' Momentum, rigidity (eV/c, T.m) : ',p, brho
      write(6,*) ' beta : ',bta
      write(6,*) ' E_tot, gamma, G.gamma : ',E, E/am, E/am * G
      write(6,*) ' beta.gamma : ', bta* E/am
      write(6,*) ' '

      write(6,*) '  Give G.gamma : '
      read(*,*,err=2) Gg

      E = Gg/G * am
      Ek = E-am
      p = sqrt(Ek*(Ek + 2.d0*am))
      brho = p/dble(nq)/c
      bta = p / sqrt(p*p + am*am)
            
      write(6,*) ' Momentum, rigidity (eV/c, T.m) : ',p, brho
      write(6,*) ' beta : ',bta
      write(6,fmt='(a,1p,3f30.8)') 
     >' E_tot, E_kin, gamma, G.gamma : ',E, Ek, E/am, E/am * G
      write(6,*) ' '

      stop
      end

 
