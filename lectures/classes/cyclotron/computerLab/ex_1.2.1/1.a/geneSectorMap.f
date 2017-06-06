      implicit double precision (a-h,o-z)
      parameter (pi = 4.d0*atan(1.d0)))

C------------ Hypothesis :
C Total angle extent of the field map
      AT = 60.d0  /180.d0*pi
C Radial extent of the field map
      Rmi = 10.d0   ! cm
      Rma = 70.d0   ! cm
C Take RM=50 cm reference radius, as this (arbitray) value is found in other exercises
      RM = 50.d0
C dR is the radial distance between two nodes, good starting point is dR = 0.5 cm
      NR = 120 + 1
      dR = (Rma - Rmi) / (NR -1)
C dX=RM*dA is the arc length between two nodes along R=RM arc, given angle increment dA
C A good starting point (by experience) is dX a few mm, say ~0.5 cm
      dX = 0.5d0   ! cm
      NX = NINT(RM*AT / dX)  +1  ! cm
      dX = RM*AT / (NX - 1)
      dA = dX / RM
      A1 = 0.d0 ; A2 = AT
C----------------------------------------------

      BX = 0.d0 ; BY = 0.d0, Z = 0.d0
      BZ = 5.d0  ! kG

      open(unit=2,file='geneSectorMap.out')
      write(2,*) '# Field map generated using geneSectorMap.f '
      write(2,fmt='(a)') '# AT,  Rmi, Rma, RM, NR, dR, NX, dX, dA '
      write(2,fmt='(a)') '# deg, cm   cm   cm      cm      cm  rd '
      write(2,fmt='(a,1p,4(e16.8,1x),2(i3,1x,e16.8,1x),e16.8)') 
     >'# ',AT, Rmi, Rma, RM, NR, dR, NX, dX, dA
      write(2,*) '# '
      write(2,*) '# R, Z, X, BY, BZ, BX '
      write(2,*) '# cm cm rd kG  kG  kG '
      write(2,*) '# '

      do jr = 1, NR
        R = Rmi + dble(NR-1)*dR
        do ix = 1, NX
          A = A1 + dble(NX-1)*dA
          write(2,fmt='(1p,6(e16.8),a)')  R, Z, A, BY, BZ, BX
        enddo
      enddo

      stop ' Job complete ! Field map is in geneSectorMap.out.'
      end
