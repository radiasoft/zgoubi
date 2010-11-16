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
C  Brookhaven National Laboratory               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE ERERZ(ER0,R,Z,IDE,
     >                             EC,DEC,DDEC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ER0(*)
      DIMENSION EC(*),DEC(2,*),DDEC(2,2,*)

      Z2 = Z*Z
      R2 = R*R
      R3 = R*R2
      R4 = R2*R2

C------------- ER0(I)=D^IER(R,Z=0)/DR^I, hence :

C----- Er(r,z)
      DERZ2 = ER0(1)/R2 - ER0(2)/R - ER0(3)
      DERZ4 = 3.D0*(-ER0(1)/R4 + ER0(2)/R3 - ER0(3)/R2)
     >  + 2.D0*ER0(4)/R + ER0(5)
      EC(1) = ER0(1) + (0.5D0 * DERZ2 + Z2/24.D0 * DERZ4)*Z2
C----- Ez(r,z)
      DEZZ = -ER0(1)/R - ER0(2)
      DEZZ3 = (ER0(1)/R - ER0(2))/R2+2.D0/R*ER0(3)+ER0(4)
      EC(2) = Z *( DEZZ + Z2/6.D0* DEZZ3 )

C-----  Derivatives
C      dEr/dr
      DECuZ2 = (-ER0(1)/R + ER0(2))/R2 - 0.5D0*(ER0(3)/R+ER0(4))
C************      DEC(1,1) = **** ER0(2)  ****  + Z2*DECuZ2
C      dEr/dz
      DEC(2,1) = Z* (DERZ2 + Z2/6.D0*DERZ4)
C      dEz/dz
      DEC(2,2) = DEZZ + 0.5D0*Z2*DEZZ3 

C      d2Er/dr2 
      DDEC(1,1,1) =.5D0*Z2*((6.D0*(ER0(1)/R2-ER0(2)/R)+3.D0*ER0(3))/R2-
     >      ER0(4)/R - ER0(5))
C      d2Er/dzdr, etc.
      DDEC(2,1,1) = 2.D0*Z*DECuZ2
      DDEC(2,2,1) = DERZ2 + 0.5D0*Z2*DERZ4                
C      DDEC(2,2,2) = **** ER0(3)  ****  + Z*DEZZ3         

      IF(IDE .EQ. 2) RETURN

      RETURN
      END
