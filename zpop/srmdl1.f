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
      SUBROUTINE SRMDL1(NLOG,KC,NPH,GPH,NPS,GPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GPH(*), GPS(*)

      CHARACTER REP

      COMMON/VXPLT/ XMI,XMA,YI,YA,KX,KY,IAX,LIS,NB
      PARAMETER (NOCT=100)

      OCT0 = GPH(1) * ( 3.D0 * (1 + GPS(1)**2) + GPH(1)**2 )/8.D0 
      OCT1 = OCT0 - GPH(1)
      OCT2 = OCT0 + GPH(1)

 2    CONTINUE

C----- Calculate min/max
      YMI = 1.D10
      YMA = -1.D10
      DO 10 IPH = 1, NPH
        PH = GPH(IPH)
        DO 11 IPS = 1, NPS
          PS = GPS(IPS)
          UPS = 1 + PS*PS
          UPS2 = UPS * UPS
          RPS = SQRT(UPS)
          HS = PH / RPS
          UU = .5D0 * HS * (3.D0 + HS * HS)
          DO 12 IOCT = 0,NOCT
            OCT = ((NOCT-IOCT)*OCT1 + IOCT*OCT2)/NOCT
            U = UU - 4.D0 * OCT / RPS**3
            ASH = ASINH(U)
            SH =  SINH ( ASINH(U) / 3.D0 )
            SH2 = SH * SH
            D = UPS2 * (1.D0 + 4.D0 * SH2 )**3
            IF(KC .EQ. 1) THEN
C------------- Sigma component
              Y = - (1.D0 - 4.D0 * SH2 ) / D
            ELSEIF(KC .EQ. 2) THEN
C------------- Pi component
              D = UPS2 * D
              Y = - 4.D0 * PS * SH / D
            ENDIF
            IF(Y .GT. YMA) THEN 
              XYM = OCT
              YMA = Y  
            ENDIF    
            IF(Y .LT. YMI) YMI = Y      
 12       CONTINUE
 11     CONTINUE
 10   CONTINUE

C      CALL TXTFBG
C      CALL LINTYP(1)
      CALL TRAXES(OCT1, OCT2, YMI,YMA,2) 

C----- Plot
      DO 20 IPH = 1, NPH
        PH = GPH(IPH)
        DO 21 IPS = 1, NPS
          PS = GPS(IPS)
          WRITE(6,*) ' Curve #', (iph-1)*nps+ips
          WRITE(6,*) '          Gamma*Psi = ',ps
          WRITE(6,*) '          Gamma*Phi = ',ph         
          UPS = 1 + PS*PS
          UPS2 = UPS * UPS
          RPS = SQRT(UPS)
          HS = PH / RPS
          UU = .5D0 * HS * (3.D0 + HS * HS)
          NPT = 0 
          DO 22 IOCT = 0,NOCT
            OCT = ((NOCT-IOCT)*OCT1 + IOCT*OCT2)/NOCT
            NPT = NPT + 1
            U = UU - 4.D0 * OCT / RPS**3
            ASH = ASINH(U)
            SH =  SINH ( ASINH(U) / 3.D0 )
            SH2 = SH * SH
            D = UPS2 * (1.D0 + 4.D0 * SH2 )**3
            IF(KC .EQ. 1) THEN
C------------- Sigma component
              Y = - (1.D0 - 4.D0 * SH2 ) / D
            ELSEIF(KC .EQ. 2) THEN
C------------- Pi component
              D = UPS2 * D
              Y = - 4.D0 * PS * SH / D
            ENDIF
            X = OCT          
            IF(OCT .EQ. OCT1) CALL VECTPL(X,Y,4)
            CALL VECTPL(X,Y,2)
            WRITE(6,*) npt, x, y
            IF(LIS .EQ. 2) CALL IMPV(NLOG,NPT,X,Y,DUM,DUM,IDUM)
 22       CONTINUE
 21     CONTINUE
 20   CONTINUE

      CALL TRKVAR(Npt,'E','(rel.)','Omgc-t','(rad)') 

      WRITE(6,*) ' X(Ymax) = ', XYM
      WRITE(6,*) ' change scales (n/y)'
      READ(5,FMT='(A1)') REP
      IF(REP .EQ. 'Y') REP = 'y'      
      IF(REP .EQ. 'y') THEN
 30     WRITE(6,*) ' GIVE XMI, XMA:'
        READ(5,*,err=30) XMI,XMA
        OCT1 = XMI
        OCT2 = XMA
 31     WRITE(6,*) ' Give Gamma*Phi:'
        READ(5,*,err=31)  GPH(1)
        GOTO 2
      ENDIF

      RETURN
      END
