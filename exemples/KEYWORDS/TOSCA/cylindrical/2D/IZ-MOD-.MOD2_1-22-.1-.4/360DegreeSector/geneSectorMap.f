      implicit double precision (a-h,o-z)
      parameter (pi = 4.d0*atan(1.d0))

C------------ Hypothesis :
C Total angle extent of the field map
      AT = 360.d0  /180.d0*pi
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
      dX = 0.5d0             ! cm mesh step at RM, approximate: allows getting NX 
      NX = NINT(RM*AT / dX)  +1  
      dX = RM*AT / DBLE(NX - 1)  ! exact mesh step at RM, corresponding to NX
      dA = dX / RM               ! corresponding delta_angle
      A1 = 0.d0 ; A2 = AT
C----------------------------------------------

      BY = 0.d0 ; BX = 0.d0 ; Z = 0.d0
      BZ = 5.d0  ! kG

      open(unit=2,file='geneSectorMap.out')
      write(2,*) Rmi,dR,dA/pi*180.d0,dZ,
     >'     !  Rmi/cm, dR/cm, dA/deg, dZ/cm'
      write(2,*) '# Field map generated using geneSectorMap.f '
      write(2,fmt='(a)') '# AT/rd,  AT/deg, Rmi/cm, Rma/cm, RM/cm,'
     >//' NR, dR/cm, NX, dX/cm, dA/rd : '
      write(2,fmt='(a,1p,5(e16.8,1x),2(i3,1x,e16.8,1x),e16.8)') 
     >'# ',AT, AT/pi*180.d0,Rmi, Rma, RM, NR, dR, NX, dX, dA
      write(2,*) '# For TOSCA:  629 121 1 22.1 1.  !IZ=1 -> 2D ; '
     >//'MOD=22 -> polar map ; .MOD2=.1 -> one map file'
      write(2,*) '# R*cosA (A:0->360), Z==0, R*sinA, BY, BZ, BX '
      write(2,*) '# cm                 cm    cm      kG  kG  kG '
      write(2,*) '# '

      do jr = 1, NR
        R = Rmi + dble(jr-1)*dR
        do ix = 1, NX
          A = A1 + dble(ix-1)*dA
C          write(2,fmt='(1p,6(e16.8),a)')  R, Z, A, BR, BZ, BA
          X = R * sin(A)
          Y = R * cos(A)
          write(2,fmt='(1p,6(e16.8),a)')  Y, Z, X, BY, BZ, BX
        enddo
      enddo

      stop ' Job complete ! Field map stored in geneSectorMap.out.'
      end
