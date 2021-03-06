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
C  Upton, NY, 11973, USA
C  -------
c     twssfc.f
c     Written after former tunesc by Frédéric Desforges
      SUBROUTINE TWSSFC(P,R, 
     >                      F0,rPARAM,C)
      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      DIMENSION P(4,4), R(6,6) 
      DIMENSION F0(6,6)
      DIMENSION C(2,2),P12(2,2),P22INV(2,2),WORK(2)
      INTEGER   IPIV,INFO
      DIMENSION IPIV(4)

C R is the local 6*6 1-turn map. P(4,4) is the corresponding transverse Edward-Teng P-matrix

      F0(1,2) = P(2,1)/P(2,2)         ! -alpha
      F0(2,1) = F0(1,2)
      F0(1,1) = P(1,1)/P(2,2)
      F0(2,2) = (P(2,1)**2+P(2,2)**2)/(P(1,1)*P(2,2))

      F0(3,4) =  P(4,3)/P(4,4)         !  -alpha
      F0(4,3) = F0(3,4) 
      F0(3,3) = P(3,3)/P(4,4)
      F0(4,4) = (P(4,3)**2+P(4,4)**2)/(P(3,3)*P(4,4))

      RPARAM  = (SQRT(P(1,1)*P(2,2))+SQRT(P(3,3)*P(4,4)))/2

      P12(1,1) = P(1,3)
      P12(1,2) = P(1,4)
      P12(2,1) = P(2,3)
      P12(2,2) = P(2,4)
      P22INV(1,1) = P(3,3)
      P22INV(1,2) = P(3,4)
      P22INV(2,1) = P(4,3)
      P22INV(2,2) = P(4,4)
      CALL DGETRF(2,2,P22INV,2,IPIV,INFO)
      CALL DGETRI(2,P22INV,2,IPIV,WORK,2,INFO)

      C = RPARAM*MATMUL(P12,P22INV)

Compute periodic dispersion from the 6*6 1-turn map
        DEN = 
     >  (-1 + R(1,1) + R(1,2)*R(2,1) + R(2,2) - R(1,1)*R(2,2) + R(1,3)
     >*R(3,1) - R(1,3)*R(2,2)*R(3,1) + R(1,2)*R(2,3)*R(3,1) + R(1,3)
     >  *R(2,1)*R(3,2) + R(2,3)*R(3,2) - 
     >    R(1,1)*R(2,3)*R(3,2) + R(3,3) - R(1,1)*R(3,3) - R(1,2)*R(2,1)
     >*R(3,3) - R(2,2)*R(3,3) + R(1,1)*R(2,2)*R(3,3) + R(1,4)*R(4,1)
     >   - R(1,4)*R(2,2)*R(4,1) + 
     >    R(1,2)*R(2,4)*R(4,1) - R(1,4)*R(2,3)*R(3,2)*R(4,1) + R(1,3)*
     >R(2,4)*R(3,2)*R(4,1) - R(1,4)*R(3,3)*R(4,1) + R(1,4)*R(2,2)*
     >  R(3,3)*R(4,1) - R(1,2)*R(2,4)*R(3,3)*R(4,1) + 
     >    R(1,3)*R(3,4)*R(4,1) - R(1,3)*R(2,2)*R(3,4)*R(4,1) + R(1,2)*
     >R(2,3)*R(3,4)*R(4,1) + R(1,4)*R(2,1)*R(4,2) + R(2,4)*R(4,2) - 
     >  R(1,1)*R(2,4)*R(4,2) + 
     >    R(1,4)*R(2,3)*R(3,1)*R(4,2) - R(1,3)*R(2,4)*R(3,1)*R(4,2) - 
     >R(1,4)*R(2,1)*R(3,3)*R(4,2) - R(2,4)*R(3,3)*R(4,2) + R(1,1)*
     >  R(2,4)*R(3,3)*R(4,2) + R(1,3)*R(2,1)*R(3,4)*R(4,2) + 
     >    R(2,3)*R(3,4)*R(4,2) - R(1,1)*R(2,3)*R(3,4)*R(4,2) + R(1,4)*
     >R(3,1)*R(4,3) - R(1,4)*R(2,2)*R(3,1)*R(4,3) + R(1,2)*R(2,4)*
     >  R(3,1)*R(4,3) + R(1,4)*R(2,1)*R(3,2)*R(4,3) + 
     >    R(2,4)*R(3,2)*R(4,3) - R(1,1)*R(2,4)*R(3,2)*R(4,3) + R(3,4)*
     >R(4,3) - R(1,1)*R(3,4)*R(4,3) - R(1,2)*R(2,1)*R(3,4)*R(4,3) -
     >   R(2,2)*R(3,4)*R(4,3) + 
     >    R(1,1)*R(2,2)*R(3,4)*R(4,3) + R(4,4) - R(1,1)*R(4,4) - R(1,2)
     >*R(2,1)*R(4,4) - R(2,2)*R(4,4) + R(1,1)*R(2,2)*R(4,4) - 
     >  R(1,3)*R(3,1)*
     >R(4,4) + R(1,3)*R(2,2)*R(3,1)*R(4,4) - 
     >    R(1,2)*R(2,3)*R(3,1)*R(4,4) - R(1,3)*R(2,1)*R(3,2)*R(4,4) - 
     >R(2,3)*R(3,2)*R(4,4) + R(1,1)*R(2,3)*R(3,2)*R(4,4) - R(3,3)*R(4,4)
     >   + R(1,1)*R(3,3)*R(4,4) + 
     >    R(1,2)*R(2,1)*R(3,3)*R(4,4) + R(2,2)*R(3,3)*R(4,4) - R(1,1)*
     >R(2,2)*R(3,3)*R(4,4))

      IF(DEN.NE. 0.D0) THEN      
        F0(1,6) = 
     > (-R(1,6) + R(1,6)*R(2,2) - R(1,2)*R(2,6) + R(1,6)*R(2,3)*R(3,2)
     > - R(1,3)*R(2,6)*R(3,2) + R(1,6)*R(3,3) - R(1,6)*R(2,2)*R(3,3) +
     >   R(1,2)*R(2,6)*R(3,3) - R(1,3)*R(3,6) + 
     >    R(1,3)*R(2,2)*R(3,6) - R(1,2)*R(2,3)*R(3,6) + R(1,6)*R(2,4)*
     >R(4,2) - R(1,4)*R(2,6)*R(4,2) - R(1,6)*R(2,4)*R(3,3)*R(4,2) + 
     >  R(1,4)*R(2,6)*R(3,3)*R(4,2) + 
     >    R(1,6)*R(2,3)*R(3,4)*R(4,2) - R(1,3)*R(2,6)*R(3,4)*R(4,2) - 
     >R(1,4)*R(2,3)*R(3,6)*R(4,2) + R(1,3)*R(2,4)*R(3,6)*R(4,2) + 
     >  R(1,6)*R(2,4)*R(3,2)*R(4,3) - 
     >    R(1,4)*R(2,6)*R(3,2)*R(4,3) + R(1,6)*R(3,4)*R(4,3) - R(1,6)*
     >R(2,2)*R(3,4)*R(4,3) + R(1,2)*R(2,6)*R(3,4)*R(4,3) - R(1,4)*
     >  R(3,6)*R(4,3) + R(1,4)*R(2,2)*R(3,6)*R(4,3) - 
     >    R(1,2)*R(2,4)*R(3,6)*R(4,3) + R(1,6)*R(4,4) - R(1,6)*R(2,2)*
     >R(4,4) + R(1,2)*R(2,6)*R(4,4) - R(1,6)*R(2,3)*R(3,2)*R(4,4) + 
     >  R(1,3)*R(2,6)*R(3,2)*R(4,4) - 
     >    R(1,6)*R(3,3)*R(4,4) + R(1,6)*R(2,2)*R(3,3)*R(4,4) - R(1,2)*
     >R(2,6)*R(3,3)*R(4,4) + R(1,3)*R(3,6)*R(4,4) - R(1,3)*R(2,2)*
     >  R(3,6)*R(4,4) + R(1,2)*R(2,3)*R(3,6)*R(4,4) - 
     >    R(1,4)*R(4,6) + R(1,4)*R(2,2)*R(4,6) - R(1,2)*R(2,4)*R(4,6) 
     >+ R(1,4)*R(2,3)*R(3,2)*R(4,6) - R(1,3)*R(2,4)*R(3,2)*R(4,6) +
     >   R(1,4)*R(3,3)*R(4,6) - 
     >    R(1,4)*R(2,2)*R(3,3)*R(4,6) + R(1,2)*R(2,4)*R(3,3)*R(4,6) -
     > R(1,3)*R(3,4)*R(4,6) + R(1,3)*R(2,2)*R(3,4)*R(4,6) - R(1,2)*
     >  R(2,3)*R(3,4)*R(4,6))/ DEN
       ELSE
        F0(1,6) = 0.D0
       ENDIF

          DEN = 
     >  (-1 + R(1,1) + R(1,2)*R(2,1) + R(2,2) - R(1,1)*R(2,2) + R(1,3)
     >*R(3,1) - R(1,3)*R(2,2)*R(3,1) + R(1,2)*R(2,3)*R(3,1) + R(1,3)*
     >  R(2,1)*R(3,2) + R(2,3)*R(3,2) - 
     >    R(1,1)*R(2,3)*R(3,2) + R(3,3) - R(1,1)*R(3,3) - R(1,2)*R(2,1)
     >*R(3,3) - R(2,2)*R(3,3) + R(1,1)*R(2,2)*R(3,3) + R(1,4)*R(4,1) 
     >  - R(1,4)*R(2,2)*R(4,1) + 
     >    R(1,2)*R(2,4)*R(4,1) - R(1,4)*R(2,3)*R(3,2)*R(4,1) + R(1,3)*
     >R(2,4)*R(3,2)*R(4,1) - R(1,4)*R(3,3)*R(4,1) + R(1,4)*R(2,2)*
     >  R(3,3)*R(4,1) - R(1,2)*R(2,4)*R(3,3)*R(4,1) + 
     >    R(1,3)*R(3,4)*R(4,1) - R(1,3)*R(2,2)*R(3,4)*R(4,1) + R(1,2)*
     >R(2,3)*R(3,4)*R(4,1) + R(1,4)*R(2,1)*R(4,2) + R(2,4)*R(4,2) -
     >   R(1,1)*R(2,4)*R(4,2) + 
     >    R(1,4)*R(2,3)*R(3,1)*R(4,2) - R(1,3)*R(2,4)*R(3,1)*R(4,2) - 
     >R(1,4)*R(2,1)*R(3,3)*R(4,2) - R(2,4)*R(3,3)*R(4,2) + R(1,1)*
     >  R(2,4)*R(3,3)*R(4,2) + R(1,3)*R(2,1)*R(3,4)*R(4,2) + 
     >    R(2,3)*R(3,4)*R(4,2) - R(1,1)*R(2,3)*R(3,4)*R(4,2) + R(1,4)*
     >R(3,1)*R(4,3) - R(1,4)*R(2,2)*R(3,1)*R(4,3) + R(1,2)*R(2,4)*
     >  R(3,1)*R(4,3) + R(1,4)*R(2,1)*R(3,2)*R(4,3) + 
     >    R(2,4)*R(3,2)*R(4,3) - R(1,1)*R(2,4)*R(3,2)*R(4,3) + R(3,4)*
     >R(4,3) - R(1,1)*R(3,4)*R(4,3) - R(1,2)*R(2,1)*R(3,4)*R(4,3) -
     >   R(2,2)*R(3,4)*R(4,3) + 
     >    R(1,1)*R(2,2)*R(3,4)*R(4,3) + R(4,4) - R(1,1)*R(4,4) - R(1,2)
     >*R(2,1)*R(4,4) - R(2,2)*R(4,4) + R(1,1)*R(2,2)*R(4,4) - R(1,3)*
     >  R(3,1)*R(4,4) + R(1,3)*R(2,2)*R(3,1)*R(4,4) - 
     >    R(1,2)*R(2,3)*R(3,1)*R(4,4) - R(1,3)*R(2,1)*R(3,2)*R(4,4) - 
     >R(2,3)*R(3,2)*R(4,4) + R(1,1)*R(2,3)*R(3,2)*R(4,4) - R(3,3)*
     >  R(4,4) + R(1,1)*R(3,3)*R(4,4) + 
     >    R(1,2)*R(2,1)*R(3,3)*R(4,4) + R(2,2)*R(3,3)*R(4,4) - R(1,1)
     >*R(2,2)*R(3,3)*R(4,4))

      IF(DEN.NE. 0.D0) THEN      
        F0(2,6) = 
     >   (-(R(1,6)*R(2,1)) - R(2,6) + R(1,1)*R(2,6) - R(1,6)*R(2,3)*
     >R(3,1) + R(1,3)*R(2,6)*R(3,1) + R(1,6)*R(2,1)*R(3,3) + R(2,6)*
     >  R(3,3) - R(1,1)*R(2,6)*R(3,3) - 
     >    R(1,3)*R(2,1)*R(3,6) - R(2,3)*R(3,6) + R(1,1)*R(2,3)*R(3,6)
     > - R(1,6)*R(2,4)*R(4,1) + R(1,4)*R(2,6)*R(4,1) + R(1,6)*R(2,4)*
     >  R(3,3)*R(4,1) - R(1,4)*R(2,6)*R(3,3)*R(4,1) - 
     >    R(1,6)*R(2,3)*R(3,4)*R(4,1) + R(1,3)*R(2,6)*R(3,4)*R(4,1) 
     >+ R(1,4)*R(2,3)*R(3,6)*R(4,1) - R(1,3)*R(2,4)*R(3,6)*R(4,1) - 
     >  R(1,6)*R(2,4)*R(3,1)*R(4,3) + 
     >    R(1,4)*R(2,6)*R(3,1)*R(4,3) + R(1,6)*R(2,1)*R(3,4)*R(4,3) 
     >+ R(2,6)*R(3,4)*R(4,3) - R(1,1)*R(2,6)*R(3,4)*R(4,3) - R(1,4)*
     >  R(2,1)*R(3,6)*R(4,3) - R(2,4)*R(3,6)*R(4,3) + 
     >    R(1,1)*R(2,4)*R(3,6)*R(4,3) + R(1,6)*R(2,1)*R(4,4) + R(2,6)
     >*R(4,4) - R(1,1)*R(2,6)*R(4,4) + R(1,6)*R(2,3)*R(3,1)*R(4,4) - 
     >  R(1,3)*R(2,6)*R(3,1)*R(4,4) - 
     >    R(1,6)*R(2,1)*R(3,3)*R(4,4) - R(2,6)*R(3,3)*R(4,4) + R(1,1)
     >*R(2,6)*R(3,3)*R(4,4) + R(1,3)*R(2,1)*R(3,6)*R(4,4) + R(2,3)*
     >  R(3,6)*R(4,4) - R(1,1)*R(2,3)*R(3,6)*R(4,4) - 
     >    R(1,4)*R(2,1)*R(4,6) - R(2,4)*R(4,6) + R(1,1)*R(2,4)*R(4,6)
     > - R(1,4)*R(2,3)*R(3,1)*R(4,6) + R(1,3)*R(2,4)*R(3,1)*R(4,6) + 
     >  R(1,4)*R(2,1)*R(3,3)*R(4,6) + 
     >    R(2,4)*R(3,3)*R(4,6) - R(1,1)*R(2,4)*R(3,3)*R(4,6) - R(1,3)
     >*R(2,1)*R(3,4)*R(4,6) - R(2,3)*R(3,4)*R(4,6) + R(1,1)*R(2,3)*
     >  R(3,4)*R(4,6))/   DEN
       ELSE
        F0(2,6) = 0.D0
       ENDIF
     
       DEN = 
     >  (-1 + R(1,1) + R(1,2)*R(2,1) + R(2,2) - R(1,1)*R(2,2) + R(1,3)
     >*R(3,1) - R(1,3)*R(2,2)*R(3,1) + R(1,2)*R(2,3)*R(3,1) + R(1,3)
     >  *R(2,1)*R(3,2) + R(2,3)*R(3,2) - 
     >    R(1,1)*R(2,3)*R(3,2) + R(3,3) - R(1,1)*R(3,3) - R(1,2)*R(2,1)
     >*R(3,3) - R(2,2)*R(3,3) + R(1,1)*R(2,2)*R(3,3) + R(1,4)*R(4,1) 
     >  - R(1,4)*R(2,2)*R(4,1) + 
     >    R(1,2)*R(2,4)*R(4,1) - R(1,4)*R(2,3)*R(3,2)*R(4,1) + R(1,3)*
     >R(2,4)*R(3,2)*R(4,1) - R(1,4)*R(3,3)*R(4,1) + R(1,4)*R(2,2)*
     >  R(3,3)*R(4,1) - R(1,2)*R(2,4)*R(3,3)*R(4,1) + 
     >    R(1,3)*R(3,4)*R(4,1) - R(1,3)*R(2,2)*R(3,4)*R(4,1) + R(1,2)*
     >R(2,3)*R(3,4)*R(4,1) + R(1,4)*R(2,1)*R(4,2) + R(2,4)*R(4,2) - 
     >  R(1,1)*R(2,4)*R(4,2) + 
     >    R(1,4)*R(2,3)*R(3,1)*R(4,2) - R(1,3)*R(2,4)*R(3,1)*R(4,2) - 
     >R(1,4)*R(2,1)*R(3,3)*R(4,2) - R(2,4)*R(3,3)*R(4,2) + R(1,1)*
     >  R(2,4)*R(3,3)*R(4,2) + R(1,3)*R(2,1)*R(3,4)*R(4,2) + 
     >    R(2,3)*R(3,4)*R(4,2) - R(1,1)*R(2,3)*R(3,4)*R(4,2) + R(1,4)*
     >R(3,1)*R(4,3) - R(1,4)*R(2,2)*R(3,1)*R(4,3) + R(1,2)*R(2,4)*
     >  R(3,1)*R(4,3) + R(1,4)*R(2,1)*R(3,2)*R(4,3) + 
     >    R(2,4)*R(3,2)*R(4,3) - R(1,1)*R(2,4)*R(3,2)*R(4,3) + R(3,4)*
     >R(4,3) - R(1,1)*R(3,4)*R(4,3) - R(1,2)*R(2,1)*R(3,4)*R(4,3) 
     >  - R(2,2)*R(3,4)*R(4,3) + 
     >    R(1,1)*R(2,2)*R(3,4)*R(4,3) + R(4,4) - R(1,1)*R(4,4) - 
     >R(1,2)*R(2,1)*R(4,4) - R(2,2)*R(4,4) + R(1,1)*R(2,2)*R(4,4) 
     >  - R(1,3)*R(3,1)*R(4,4) + R(1,3)*R(2,2)*R(3,1)*R(4,4) - 
     >    R(1,2)*R(2,3)*R(3,1)*R(4,4) - R(1,3)*R(2,1)*R(3,2)*R(4,4) 
     >- R(2,3)*R(3,2)*R(4,4) + R(1,1)*R(2,3)*R(3,2)*R(4,4) - 
     >  R(3,3)*R(4,4) + R(1,1)*R(3,3)*R(4,4) + 
     >    R(1,2)*R(2,1)*R(3,3)*R(4,4) + R(2,2)*R(3,3)*R(4,4) - R(1,1)
     >*R(2,2)*R(3,3)*R(4,4))

      IF(DEN.NE. 0.D0) THEN      
             F0(3,6) = 
     >   (-(R(1,6)*R(3,1)) + R(1,6)*R(2,2)*R(3,1) - R(1,2)*R(2,6)*R(3,1)
     > - R(1,6)*R(2,1)*R(3,2) - R(2,6)*R(3,2) + R(1,1)*R(2,6)*R(3,2)
     >   - R(3,6) + R(1,1)*R(3,6) + 
     >    R(1,2)*R(2,1)*R(3,6) + r(2,2)*r(3,6) - r(1,1)*r(2,2)*r(3,6) -
     > r(1,6)*r(2,4)*r(3,2)*r(4,1) + r(1,4)*r(2,6)*r(3,2)*r(4,1) - 
     >  r(1,6)*r(3,4)*r(4,1) + 
     >    r(1,6)*r(2,2)*r(3,4)*r(4,1) - r(1,2)*r(2,6)*r(3,4)*r(4,1) + 
     >r(1,4)*r(3,6)*r(4,1) - r(1,4)*r(2,2)*r(3,6)*r(4,1) + r(1,2)*
     >  r(2,4)*r(3,6)*r(4,1) + r(1,6)*r(2,4)*r(3,1)*r(4,2) - 
     >    r(1,4)*r(2,6)*r(3,1)*r(4,2) - r(1,6)*r(2,1)*r(3,4)*r(4,2) - 
     >r(2,6)*r(3,4)*r(4,2) + r(1,1)*r(2,6)*r(3,4)*r(4,2) + r(1,4)*
     >  R(2,1)*R(3,6)*R(4,2) + R(2,4)*R(3,6)*R(4,2) - 
     >    R(1,1)*R(2,4)*R(3,6)*R(4,2) + R(1,6)*R(3,1)*R(4,4) - R(1,6)*
     >R(2,2)*R(3,1)*R(4,4) + R(1,2)*R(2,6)*R(3,1)*R(4,4) + R(1,6)*
     >  R(2,1)*R(3,2)*R(4,4) + R(2,6)*R(3,2)*R(4,4) - 
     >    R(1,1)*R(2,6)*R(3,2)*R(4,4) + R(3,6)*R(4,4) - R(1,1)*R(3,6)*
     >R(4,4) - R(1,2)*R(2,1)*R(3,6)*R(4,4) - R(2,2)*R(3,6)*R(4,4) +
     >   R(1,1)*R(2,2)*R(3,6)*R(4,4) - 
     >    R(1,4)*R(3,1)*R(4,6) + R(1,4)*R(2,2)*R(3,1)*R(4,6) - R(1,2)*
     >R(2,4)*R(3,1)*R(4,6) - R(1,4)*R(2,1)*R(3,2)*R(4,6) - R(2,4)*
     >  R(3,2)*R(4,6) + R(1,1)*R(2,4)*R(3,2)*R(4,6) - 
     >    R(3,4)*R(4,6) + R(1,1)*R(3,4)*R(4,6) + R(1,2)*R(2,1)*R(3,4)*
     >R(4,6) + R(2,2)*R(3,4)*R(4,6) - R(1,1)*R(2,2)*R(3,4)*R(4,6)) 
     >  /   DEN
      ELSE
             F0(3,6) = 0.D0
      ENDIF
       
      DEN =  
     >  (((R(3,2)*R(4,1) - R(3,1)*R(4,2))*(R(2,3)*R(4,1) - R(2,1)*
     >R(4,3)) - (-R(4,1) + R(2,2)*R(4,1) - R(2,1)*R(4,2))*(-R(4,1) + 
     >  R(3,3)*R(4,1) - R(3,1)*R(4,3)))*
     >     ((R(3,2)*R(4,1) - R(3,1)*R(4,2))*(-1 + R(1,1) + R(1,4)*
     >R(4,1) + R(4,4) - R(1,1)*R(4,4)) - 
     >       (R(1,2)*R(4,1) + R(4,2) - R(1,1)*R(4,2))*(R(3,1) + 
     >R(3,4)*R(4,1) - R(3,1)*R(4,4))) - 
     >    ((R(3,2)*R(4,1) - R(3,1)*R(4,2))*(R(1,3)*R(4,1) + R(4,3)
     > - R(1,1)*R(4,3)) - (R(1,2)*R(4,1) + R(4,2) - R(1,1)*R(4,2))*(
     >  -R(4,1) + R(3,3)*R(4,1) - R(3,1)*R(4,3)))*
     >     ((R(3,2)*R(4,1) - R(3,1)*R(4,2))*(R(2,1) + R(2,4)*R(4,1)
     > - R(2,1)*R(4,4)) - (-R(4,1) + R(2,2)*R(4,1) - R(2,1)*R(4,2))*
     >  (R(3,1) + R(3,4)*R(4,1) - R(3,1)*R(4,4))))

      IF(DEN.NE. 0.D0) THEN      
             F0(4,6) = 
     >      (((R(3,2)*R(4,1) - R(3,1)*R(4,2))*(R(2,3)*R(4,1) - R(2,1)*
     >R(4,3)) - (-R(4,1) + R(2,2)*R(4,1) - R(2,1)*R(4,2))*(-R(4,1) + 
     >  R(3,3)*R(4,1) - R(3,1)*R(4,3)))*
     >     ((R(3,2)*R(4,1) - R(3,1)*R(4,2))*(-(R(1,6)*R(4,1)) - R(4,6)
     > + R(1,1)*R(4,6)) - (R(1,2)*R(4,1) + R(4,2) - R(1,1)*R(4,2))*(
     >  -(R(3,6)*R(4,1)) + R(3,1)*R(4,6))) -
     >  ((R(3,2)*R(4,1) - R(3,1)*R(4,2))*(R(1,3)*R(4,1) + R(4,3) 
     >- R(1,1)*R(4,3)) - 
     >       (R(1,2)*R(4,1) + R(4,2) - R(1,1)*R(4,2))*(-R(4,1) + R(3,3)
     >*R(4,1) - R(3,1)*R(4,3)))*
     >     ((R(3,2)*R(4,1) - R(3,1)*R(4,2))*(-(R(2,6)*R(4,1)) + R(2,1)
     >*R(4,6)) - (-R(4,1) + R(2,2)*R(4,1) - R(2,1)*R(4,2))*(-(R(3,6)*
     >  R(4,1)) + R(3,1)*R(4,6))))/   DEN
      ELSE
             F0(4,6) = 0.D0 
      ENDIF

      RETURN                 

      END
