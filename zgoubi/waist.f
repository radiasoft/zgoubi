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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE WAIST
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      COMMON/INIT/ FA0(6,6),FA1(6,6),BID(6),BID1(6),IF
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
 
      DIMENSION F2(MXT),F3(MXT)
      INTEGER HV
 
 
      JDMAX=IDMAX
      JMAXT=IMAXT
 
      HV = 0
 10   CONTINUE
      HV = HV+1
 
      IF    ( HV .EQ. 1 ) THEN
        DO 4 I=1,IMAX
          IF( IEX(I) .LT. -1 ) GOTO 4
          F2(I)=F(2,I)
          F3(I)=F(3,I)
 4      CONTINUE
      ELSEIF( HV .EQ. 2 ) THEN
        DO 5 I=1,IMAX
          IF( IEX(I) .LT. -1) GOTO 5
          F2(I)=F(4,I)
          F3(I)=F(5,I)
 5      CONTINUE
      ENDIF
 
      SWI2=0D0
      SWTI2=0D0
      DO 3 ID=1,JDMAX
        IMAX1=1+(ID-1)*JMAXT
        IMAX2=IMAX1+JMAXT-1
 
        IMAXI = 0
        YMI=0D0
        YMIN= 1.D10
        YMAX=-1.D10
        WI=0D0
        YIO=F2(IMAX1)
        TMI=0D0
        TMIN= 1.D10
        TMAX=-1.D10
        WTI=0D0
        TIO=F3(IMAX1)
        DO 2 I=IMAX1,IMAX2
          IF( IEX(I) .LT. -1) GOTO 2
          IMAXI = IMAXI+1
 
          YI=F2(I) -YIO
          YMI=YMI+YI
          IF(YI .LT. YMIN) YMIN=YI
          IF(YI .GT. YMAX) YMAX=YI
          WI=WI+YI*YI
 
          TI=F3(I) -TIO
          TMI=TMI+TI
          IF(TI .LT. TMIN) TMIN=TI
          IF(TI .GT. TMAX) TMAX=TI
          WTI=WTI+TI*TI
 2      CONTINUE
 
        YMI=YMI/(DBLE(IMAXI))
        WI=   ( WI/DBLE(IMAXI) - YMI*YMI )*5.5225D0
        SWI2 = SWI2 + WI
        WI = SQRT(WI)
        WT = YMAX-YMIN
 
        TMI=TMI/(DBLE(IMAXI))
        WTI=   ( WTI/DBLE(IMAXI) - TMI*TMI )*5.5225D0
        SWTI2 = SWTI2 + WTI
        WTI = SQRT(WTI)
        WTT = TMAX-TMIN
 
    3 CONTINUE
 
      IF    ( HV .EQ. 1 ) THEN
        FA1(1,1) = SQRT(SWI2 /JDMAX)
        FA1(2,2) = SQRT(SWTI2/JDMAX) 
        FA1(2,1) = SQRT(FA1(1,1)*FA1(2,2)-1.D0)
C------- sign -alpha
C        IF() FA1(2,1) =-FA1(2,1) 
        FA1(1,2) =FA1(2,1) 
      ELSEIF( HV .EQ. 2 ) THEN
        FA1(3,3) = SQRT(SWI2 /JDMAX)
        FA1(4,4) = SQRT(SWTI2/JDMAX) 
        FA1(3,4) = SQRT(FA1(3,3)*FA1(4,4)-1.D0)
C------- sign -alpha
C        IF() FA1(4,3) =-FA1(4,3) 
        FA1(3,4) =FA1(4,3) 
      ENDIF
 
      IF( HV .EQ. 1 ) GOTO 10
 
      RETURN
      END
