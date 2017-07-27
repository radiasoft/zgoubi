      implicit double precision (a-h,o-z)
      parameter (c=2.99792458d8) 

      am = 238.0289 * 931.494061d6
      nq = 92
      nuc = 238
      write(6,*) ' U mass, #q, #nuc : ', am , nq, nuc
      write(6,*) ' q/A = ',nq,'/',nuc,'  = ',dble(nq) /dble(nuc) 

 2    continue
      write(6,*) '  Give kinetic E per nucleon (eV) : '
      read(*,*,err=2) Ekn
      Ek = Ekn * dble(nuc)

      p = sqrt(Ek*(Ek + 2.d0*am))
      brho = p/dble(nq)/c
      bta = p / sqrt(p*p + am*am)
            
      write(6,*) ' Momentum, rigidity (eV/c, T.m) : ',p, brho
      write(6,*) ' beta : ',bta
      write(6,*) ' '

      stop
      end

 
