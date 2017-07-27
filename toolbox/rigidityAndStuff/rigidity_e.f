      implicit double precision (a-h,o-z)
      parameter (c=2.99792458d8) 

C________________________________

      am = 0.51099892d6
      G = 1.15965218076d-3
      write(6,*) '  electron '

         goto 12

 14      continue
      write(6,*) '  give electron momentum (eV/c) : '
      read(*,*,err=15) p

      bta = p / sqrt(p*p -am*am)
      brho = p / c
      E = sqrt(p*p + am*am)
      Ek = E - am
            
      write(6,*) ' Energie, total, kinetic (eV) : ',E,Ek
      write(6,*) ' Rigidity (T.m) : ',Brho
      write(6,*) ' Normalized : rigidity/33.35640843 : ',
     >      Brho/33.35640843
      write(6,*) ' beta : ',bta
      write(6,*) ' G.gamma : ',G*E/am
      write(6,*) ' '
      write(6,*) ' '

        goto 14

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
      write(6,*) ' G.gamma : ',G*E/am
      write(6,*) ' '
      write(6,*) ' '

        goto 15

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
      write(6,*) ' G.gamma : ',G*E/am
      write(6,*) ' '
      write(6,*) ' '

 12   continue
      write(6,*) '  give total E (eV)  : '
      read(*,*,err=12) E

      Ek = E - am
      p = sqrt(Ek*(Ek + 2.d0*am))
      brho = p/c
      bta = p / sqrt(p*p + am*am)
            
      write(6,*) ' momentum, rigidity (eV/c, T.m) : ',p, brho
      write(6,*) ' Brho / Brho[E+M=17.198d9 eV] : ', brho/57.366353
      write(6,*) ' G.gamma : ',G*E/am
      write(6,*) ' beta : ',bta

        goto 12

      stop
      end

 
