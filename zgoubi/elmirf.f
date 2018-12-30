C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  François Méot
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.
C
C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 51 Franklin Street, Fifth Floor,
C  Boston, MA  02110-1301  USA
C
C  François Méot <fmeot@bnl.gov>
C  Brookhaven National Laboratory
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE ELMIRF(XX,Z,BRI,
     >                           E,DE,DDE)
C Mirror, straight slits. Frame is Cartesian
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(5,3),DE(3,3),DDE(3,3,3)

      INCLUDE "C.AIM.H"     ! COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "C.ORDRES.H"     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ

C----- D=gap   AL=plate lengths, V=potentials
      DIMENSION AL(7), V(7), ALIN(7), VIN(7)
      SAVE D, NP, AL, V

      DIMENSION DV(7),
     >      CHDX(7),CHDX2(7),CHDX3(7),SHDX(7),SHDX2(7),SHDX3(7)

      DO 20 I=2, NP
 20     DV(I) = (V(I)-V(I-1))*BRI

C------- 3-plate electrical mirror potential. For checks :
C      VR =
C     >  V21D*D*ATAN(SINH(PI* X    /D))/COS(PI*Z/D) +
C     >  V32D*D*ATAN(SINH(PI*(X-AL)/D))/COS(PI*Z/D)
C      CALL ENRGW(X,VR)
C--------------------------------

      PID = PI / D
      PID2 = PID*PID
      SECDZ = 1.D0 / COS(PID*Z)
      SECDZ2 = SECDZ*SECDZ
      SECDZ3 =SECDZ2*SECDZ
      SECDZ5 = SECDZ3*SECDZ2
      TDZ = TAN(PID*Z)

      X = XX
      DO 21 I=2, NP
        X = X - AL(I-1)
        CHDX(I)= DCOSH(PID*(X+AL(1)))
        CHDX2(I) = CHDX(I)*CHDX(I)
        CHDX3(I) =CHDX2(I)*CHDX(I)
        SHDX (I)= DSINH(PID*(X+AL(1)))
        SHDX2(I) = SHDX(I)*SHDX(I)
 21     SHDX3(I) = SHDX2(I)*SHDX(I)

C----- EX = - DV/DX
      E(1,1)=0.D0
      DO 22 I=2, NP
 22      E(1,1)= E(1,1) +
     >     ( - (DV(I)*CHDX(I)*SECDZ)/(1.D0 + SECDZ2*SHDX2(I)) )

C----- EZ = - DV/DZ
      E(1,3)=0.D0
      DO 23 I=2, NP
 23     E(1,3)= E(1,3) +
     >  (- (DV(I)*SECDZ*SHDX(I)*TDZ)/ (1.D0 + SECDZ2*SHDX2(I)) )

C     ... DEX/DX
      DE(1,1) =0.D0
      DO 24 I=2, NP
 24     DE(1,1) = DE(1,1)+
     >   ((2.D0*PID*DV(I)*CHDX2(I)*SECDZ3*SHDX(I))
     >                                       /(1.D0+SECDZ2*SHDX2(I))**2-
     >     (PID*DV(I)*SECDZ*SHDX(I)) /(1.D0 + SECDZ2*SHDX2(I)) )

C     .. DEX/DZ = DEZ/DX
      DE(3,1) = 0.D0
      DO 25 I=2, NP
 25     DE(3,1) = DE(3,1) +
     >  ( (2.D0*PID*DV(I)*CHDX(I)*SECDZ3*SHDX2(I)*TDZ)
     >                                   /(1.D0+SECDZ2*SHDX2(I))**2 -
     >   (PID*DV(I)*CHDX(I)*SECDZ*TDZ)/(1.D0 + SECDZ2*SHDX2(I)) )

C     ... D2EX/DX2
      DDE(1,1,1)= 0.D0
      DO 26 I=2, NP
 26      DDE(1,1,1)= DDE(1,1,1)+
     >   (- (8.D0*PID2*DV(I)*CHDX3(I)*SECDZ5*SHDX2(I))
     >                                     /(1.D0+SECDZ2*SHDX2(I))**3+
     >   (2.D0*PID2*DV(I)*CHDX3(I)*SECDZ3)/(1.D0 + SECDZ2*SHDX2(I))**2 +
     >   (6.D0*PID2*DV(I)*CHDX(I)*SECDZ3*SHDX2(I))
     >                                    /(1.D0+SECDZ2*SHDX2(I))**2-
     >   (PID2*DV(I)*CHDX(I)*SECDZ)/(1.D0 + SECDZ2*SHDX2(I)) )

C     .. D2EX/DXDZ = D2EZ/DX2
      DDE(3,1,1)=  0.D0
      DO 27 I=2, NP
 27      DDE(3,1,1)=  DDE(3,1,1) +
     >  (-(8.D0*PID2*DV(I)*CHDX2(I)*SECDZ5*SHDX3(I)*TDZ)
     >                                 /(1.D0 + SECDZ2*SHDX2(I))**3 +
     >   (6.D0*PID2*DV(I)*CHDX2(I)*SECDZ3*SHDX(I)*TDZ)
     >                                  /(1.D0 + SECDZ2*SHDX2(I))**2 +
     >   (2.D0*PID2*DV(I)*SECDZ3*SHDX3(I)*TDZ)
     >                                 /(1.D0 + SECDZ2*SHDX2(I))**2 -
     >   (PID2*DV(I)*SECDZ*SHDX(I)*TDZ)/(1.D0 + SECDZ2*SHDX2(I)) )

      RETURN

      ENTRY ELMIRW(DIN,ALIN,VIN,NPIN)
      D = DIN
      NP=NPIN
      DO 10 I=1,NP
        AL(I) = ALIN(I)
 10     V(I) = VIN(I)

      RETURN
      END
