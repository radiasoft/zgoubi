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
C  USA
C  -------
      SUBROUTINE INIDRT(TITL,ND,XI,
     >                             RCS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) TITL
      DIMENSION RCS(*)

      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.DROITE.H"     ! COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT

      INTEGER DEBSTR

      IDRT =A(NOEL,ND-20) 
 
      IF    (IDRT .EQ. 1) THEN
        CA(2) =  A(NOEL,31)
        SA(2) =  A(NOEL,32)
        CM(2) =  A(NOEL,33)
        RCS(2) = SQRT(CA(2)*CA(2)+SA(2)*SA(2))
        CA(2) =  CA(2)/RCS(2)
        SA(2) =  SA(2)/RCS(2)
        CM(2) =  CM(2)/RCS(2)
      ELSEIF(IDRT .EQ. -1) THEN
        CA(1) =  A(NOEL,31)
        SA(1) =  A(NOEL,32)
        CM(1) =  A(NOEL,33)
        RCS(1) = SQRT(CA(1)*CA(1)+SA(1)*SA(1))
        CA(1) =  CA(1)/RCS(1)
        SA(1) =  SA(1)/RCS(1)
        CM(1) =  CM(1)/RCS(1)
      ELSEIF(IDRT .GE. 2) THEN
        II = -3
        DO 1 I = 1, IDRT
          II = II + 3
          CA(I) =  A(NOEL,30+II+1)
          SA(I) =  A(NOEL,30+II+2)
          CM(I) =  A(NOEL,30+II+3)
          RCS(I) = SQRT(CA(I)*CA(I)+SA(I)*SA(I))
          CA(I) =  CA(I)/RCS(I)
          SA(I) =  SA(I)/RCS(I)
          CM(I) =  CM(I)/RCS(I)
 1      CONTINUE 
      ENDIF

         IF(TITL(DEBSTR(TITL):DEBSTR(TITL)+5) .EQ. 'MIRROR')  THEN
C---------- Allow backward ray-tracing in maps of mirror lenses
C----------
           CALL TRANSW(.TRUE.,.TRUE.)
           IF    (IDRT.EQ.0) THEN
             CA(1)=1.D0
             SA(1)=0.D0
             CM(1)= -XI
             RCS(1) = 1.D0
             CA(2)=1.D0
             SA(2)=0.D0
             CM(2)= -XI
             RCS(2) = 1.D0
             IDRT = 2
           ELSEIF(IDRT.EQ.-1) THEN
             CA(2)=1.D0
             SA(2)=0.D0
             CM(2)= -XI
             RCS(2) = 1.D0
             IDRT = 2
           ELSEIF(IDRT.EQ.1) THEN
             CA(1)=1.D0
             SA(1)=0.D0
             CM(1)= -XI
             RCS(1) = 1.D0
             IDRT = 2
           ENDIF
         ELSEIF(TITL(DEBSTR(TITL):DEBSTR(TITL)+5) .EQ. 'BACK')  THEN
C---------- Allow backward ray-tracing in maps such as SPES3
           CALL TRANSW(.FALSE.,.TRUE.)
         ENDIF
        
      RETURN
      END
