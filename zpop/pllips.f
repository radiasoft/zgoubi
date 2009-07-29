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
      SUBROUTINE PLLIPS(NLOG,LM,KX,KY,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH
      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7
      COMMON/INPVR/ KVAR, KPOL, KDIM
      INCLUDE 'MAXNTR.H'
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR

      DIMENSION YM(3), YPM(3), U(3), A(3), B(3)
      DIMENSION YMX(6), YPMX(6)
      DIMENSION XSIGU(3)
      CHARACTER TXT*80

      IF(NPTR .GT. 0) THEN
        IF(NPTS.GT. NPTR) NPTS=NPTR

      KWRIT=0
      CALL LPSFIT(NLOG,KWRIT,LM,
     >                          YM,YPM,YMX,YPMX,U,A,B,*60,*60)
 60   CONTINUE

      DO 1 I=1,3
 1      XSIGU(I)=1.D0
      CALL LPSCNT(YM,YPM,U,A,B,XSIGU,NLOG,
     >                                    NCOUNT)
      IF    (KX.EQ.2 .OR. KX.EQ.12) THEN
        JJ=1
      ELSEIF(KX.EQ.4 .OR. KX.EQ.14) THEN
        JJ=2
      ELSE
        JJ=3
      ENDIF

          CALL LINTYP(1)
          R= SQRT(XSIGU(JJ) * U(JJ) * B(JJ)) 
          Y = R + YM(JJ)
          YP= -A(JJ)/B(JJ) * R + YPM(JJ)
          CALL VECTPL(Y,YP,4)      
          DO 201 IPHI=0,400
            PHI = IPHI*DPI/400
            Y = R * COS(PHI)
            YP= (-A(JJ)* Y + R * SIN(PHI) )/B(JJ)  + YPM(JJ)
            Y = Y + YM(JJ)
            CALL VECTPL(Y,YP,2)      
            IF(LIS.EQ.2) CALL IMPV(NLOG,0,Y,YP,DUM,DUM,IDUM,IDUM,KX,KY)
 201      CONTINUE

          CALL FBGTXT
          WRITE(*,FMT='(/,'' MATCHING ELLIPSE : '')')
          WRITE(*,108) KVAR(KY),KVAR(KX)
 108      FORMAT(A,' vs. ',A) 
          J = 2*JJ - 1
          WRITE(*,107) YMX(J),YMX(J+1),YPMX(J),YPMX(J+1)
 107      FORMAT(' Min-max - Hor.: ',1P,2G13.5,'; Ver.: ',2G13.5)
          WRITE(*,106) YM(J),YPM(J)
 106      FORMAT(' Centering,  H/V : ',1P,G11.3,' / ',G11.3)
          WRITE(*,105) U(JJ),B(JJ),A(JJ)
 105      FORMAT(' Eps/pi, Beta, Alpha: ',1P,G12.4,2G11.3)

          WRITE(TXT,105) U(JJ),B(JJ),A(JJ)
          CALL TRTXT(10.D0,18.D0,TXT,80,0)
          CALL LINTYP(-1)

      ENDIF
      RETURN 
C      RETURN 1
      END
