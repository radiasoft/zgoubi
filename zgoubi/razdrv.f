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
      SUBROUTINE RAZDRV(IOP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     --------------------------------------------------------
C     INITIALISATIONS DANS CHAMC.
C     INITIALISATIONS LIEES A L'ORDRE DE CALCUL ( 2, 25 OU 4 )
C     DU Champ B(X,Y,Z).
C     --------------------------------------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CHAMP.H"     ! COMMON/CHAMP/ BZ0(5,5), EZ0(5,5)
      INCLUDE "C.CHAVE.H"     ! COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      INCLUDE "C.DDBXYZ_2.H"     ! COMMON/DDBXYZ/ DB(9),DDB(27)
      INCLUDE "C.D3B.H"     ! COMMON/D3BXYZ/ D3BX(27), D3BY(27), D3BZ(27)
      INCLUDE "C.D4B_2.H"     ! COMMON/D4BXYZ/ D4BX(81) ,D4BY(81) ,D4BZ(81)
      INCLUDE "C.DDEXYZ_2.H"     ! COMMON/DDEXYZ/ DE(9),DDE(27)
      INCLUDE "C.D3E.H"     ! COMMON/D3EXYZ/ D3EX(27), D3EY(27), D3EZ(27)
      INCLUDE "C.D4EXYZ_2.H"     ! COMMON/D4EXYZ/ D4EX(81), D4EY(81), D4EZ(81)

      GOTO (1,2,1) IOP
 
 1    CONTINUE

      DO 10 I=1,3
 10     B(1,I)=0.D0
        DO 11 I=1,9
 11       DB(I)=0.D0
          DO 12 I=1,27
 12         DDB(I)=0.D0
            DO 13 I=1,27
              D3BX(I)=0.D0
              D3BY(I)=0.D0
              D3BZ(I)=0.D0
 13         CONTINUE
              DO 14 I=1,81
                D4BX(I)=0.D0
                D4BY(I)=0.D0
                D4BZ(I)=0.D0
 14           CONTINUE

      DO 15 J=1,5
        DO 15 I=1,5
 15       BZ0(I,J) = 0.D0

      IF(IOP .EQ. 1) RETURN
 
 2    CONTINUE

      DO 20 I=1,3
 20     E(1,I)=0.D0
        DO 21 I=1,9
 21       DE(I)=0.D0
          DO 22 I=1,27
 22         DDE(I)=0.D0
            DO 23 I=1,27
              D3EX(I)=0.D0
              D3EY(I)=0.D0
              D3EZ(I)=0.D0
 23         CONTINUE
              DO 24 I=1,81
                D4EX(I)=0.D0
                D4EY(I)=0.D0
                D4EZ(I)=0.D0
 24           CONTINUE

      DO 25 J=1,5
        DO 25 I=1,5
 25       EZ0(I,J) = 0.D0

      RETURN
      END
