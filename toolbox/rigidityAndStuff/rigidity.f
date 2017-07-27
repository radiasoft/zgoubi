      implicit double precision (a-h,o-z)
      parameter (c=2.99792458d8) 

      am = 938.27231d6
      write(6,*) '  proton mass and q/=1 assumed '

  4   continue
      write(6,*) '  give G.gamma : '
      read(*,*,err=14) Gg

      gma = Gg / 1.7928474
      E = gma * am
      p = sqrt(E*E - am*am)
      brho = p / c
      Ek = E - am
      bta = p / sqrt(p*p + am*am)
            
      write(6,*) ' Energie, total, kinetic (eV) : ',E,Ek
      write(6,*) ' Rigidity (T.m) : ',Brho
      write(6,*) ' beta : ',bta
      write(6,*) ' '
      write(6,*) ' '
      write(6,*) ' '

 14   continue
      write(6,*) '  give beta proton : '
      read(*,*,err=14) bta

      p = am*bta / sqrt(1.d0-bta*bta)
      brho = p / c
      E = sqrt(p*p + am*am)
      Ek = E - am
            
      write(6,*) ' Energie, total, kinetic (eV) : ',E,Ek
      write(6,*) ' Rigidity (T.m) : ',Brho
      write(6,*) ' beta : ',bta
      write(6,*) ' '
      write(6,*) ' '
      write(6,*) ' '

 1    continue
      write(6,*) '  give momentum (eV/c) (q/=1 assumed) : '
      read(*,*,err=1) p

      brho = p / c
      E = sqrt(p*p + am*am)
      Ek = E - am
      bta = p / sqrt(p*p + am*am)
            
      write(6,*) ' Energie, total, kinetic (eV) : ',E,Ek
      write(6,*) ' Rigidity (T.m) : ',Brho
      write(6,*) ' beta : ',bta
      write(6,*) ' '
      write(6,*) ' '
      write(6,*) ' '

 2    continue
      write(6,*) '  give kinetic E, eV (q/=1 assumed) : '
      read(*,*,err=2) Ek

      p = sqrt(Ek*(Ek + 2.d0*am))
      brho = p/c
      bta = p / sqrt(p*p + am*am)
            
      write(6,*) ' momentum, rigidity (eV/c, T.m) : ',p, brho
      write(6,*) ' beta : ',bta
      write(6,*) ' '

 3    continue
      write(6,*) '  give kinetic E, eV (q/=1 assumed) : '
      read(*,*,err=2) Ek

      p = sqrt(Ek*(Ek + 2.d0*am))
      brho = p/c
      bta = p / sqrt(p*p + am*am)
            
      write(6,*) ' momentum, rigidity (eV/c, T.m) : ',p, brho
      write(6,*) ' beta : ',bta
      write(6,*) ' '

C________________________________

      am = 0.511e6
      write(6,*) '  electron '

 15   continue
      write(6,*) '  give beta electron : '
      read(*,*,err=15) bta

      p = am*bta / sqrt(1.d0-bta*bta)
      brho = p / c
      E = sqrt(p*p + am*am)
      Ek = E - am
            
      write(6,*) ' Energie, total, kinetic (eV) : ',E,Ek
      write(6,*) ' Rigidity (T.m) : ',Brho
      write(6,*) ' beta : ',bta
      write(6,*) ' '
      write(6,*) ' '
      write(6,*) ' '

 11   continue
      write(6,*) '  give rigidity, T.m (q/=1 assumed) : '
      read(*,*,err=11) brho

      p = brho * c
      E = sqrt(p*p + am*am)
      Ek = E - am
      bta = p / sqrt(p*p + am*am)
            
      write(6,*) ' Energie, total, kinetic (eV) : ',E,Ek
      write(6,*) ' Rigidity (T.m) : ',Brho
      write(6,*) ' beta : ',bta
      write(6,*) ' '
      write(6,*) ' '

 12   continue
      write(6,*) '  give kinetic E, eV (q/=1 assumed) : '
      read(*,*,err=12) Ek

      p = sqrt(Ek*(Ek + 2.d0*am))
      brho = p/c
            
      write(6,*) ' momentum, rigidity (eV/c, T.m) : ',p, brho
      stop
      end

 
