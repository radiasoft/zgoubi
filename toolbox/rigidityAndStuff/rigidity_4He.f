      implicit double precision (a-h,o-z)
      parameter (c=2.99792458d8) 

      am = 2808.39148d6  /3.*4     ! approximated from 3He, make it exact 
      Q = 2.
      A = 4  ! nucleons
      write(6,*) '  Mass, Q, Q/A, A : ',am, Q, Q/A, A

 14   continue
      write(6,*) '  Give rigidity : '
      read(*,*,err=14) brho

      p = Q * brho * c
      E = sqrt(p*p + am*am)
      Ek = E - am
            
      write(6,*) 'Rigidity : ',brho
      write(6,*) 'Energie, total, kinetic (eV/A) : ',E/A,Ek/A
      write(6,*) 'beta, beta*gamma, Ggamma : ',bta, bta*E/am, G*E/am
      write(6,*) ' '
      write(6,*) ' '
      write(6,*) ' '

      write(6,*) '  Give kinetic E (eV/nucleon) : '
      read(*,*,err=14) EkA
      Ek = EkA * A
      p = sqrt(Ek*(Ek+2.*am))
      brho = p / (Q*c)
      E = Ek + am
            
      write(6,*) 'Rigidity : ',brho
      write(6,*) 'Energie/N, total, kinetic (eV/A) : ',E/A,Ek/A
      write(6,*) 'beta, beta*gamma : ',bta, bta*E/am
      write(6,*) ' '
      write(6,*) ' '
      write(6,*) ' '

      stop
      end

 
