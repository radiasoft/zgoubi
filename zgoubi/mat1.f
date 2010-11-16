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
      SUBROUTINE MAT1(R,T,IT1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*) , T(6,6,*)
C     -------------------------------------------------
C     OPTION  IORD = 1 :
C       MATRICES ORDRE 1 ET 2, SERIES DE TAYLOR ORDRE 3
C     -------------------------------------------------
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
 
      CALL RAZ(R,6*6)
      R(5,5) = 1.D0
      R(6,6) = 1.D0
      CALL RAZ(T,5*6*6)
      S1 = 2.D0 * F(6,IT1)
 
C..............................................
C             Y     T     Z     P     L     D
C
C        Y   R11   R12   R13   R14   R15   R16
C        T   R21   R22   R23   R24   R25   R26
C        Z   R31   R32   R33   R34   R35   R36
C        P   R41   R42   R43   R44   R45   R46
C        L   R51   R52   R53   R54   R55   R56
C        D   R61   R62   R63   R64   R65   R66
C..............................................

      I10 = IT1+9
      I11 = IT1+10
      DP = ( FO(1,I10) - FO(1,I11) ) / (.5D0*( FO(1,I10) + FO(1,I11) ) )
      DP2 = DP*DP
      DO 11 J=2,5
        R(J-1,6)  = (F(J,I10) - F(J,I11)) /DP
        DO 11 I=1,4
          I2 = 2*I +IT1-1
          I3 = I2+1 
          UO = FO(I+1,I2)-FO(I+1,I3)
          R(J-1,I)  = (F(J,I2) - F(J,I3)) / UO
          IF(J .EQ. 5) THEN
            R(5,I)  = ( F(6,I2) - F(6,I3) ) / UO
          ENDIF
C          write(*,*) j,I2,I3, F(J,I2), F(J,I3)
 11   CONTINUE
      R(5,6)  = ( F(6,I10) - F(6,I11) ) / DP
 
      IF(IMAX.NE.13) RETURN
C     ... Compute Ri5, i=1,5. 
      IF(FO(6,13)-FO(6,12) .NE. 0.D0) THEN
        DL=FO(6,13)-FO(6,12)
        R(1,5)  = ( F(2,13) - F(2,12) ) / DL
        R(2,5)  = ( F(3,13) - F(3,12) ) / DL
        R(3,5)  = ( F(4,13) - F(4,12) ) / DL
        R(4,5)  = ( F(5,13) - F(5,12) ) / DL
        R(5,5)  = ( F(6,13) - F(6,12) ) / DL
      ENDIF
      RETURN
      END
