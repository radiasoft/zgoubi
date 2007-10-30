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
      SUBROUTINE CLASS(NRM,NC,CX,CY,
     >                              XFB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CX(*), CY(*)

      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      PARAMETER(MXC=400)      
      DIMENSION CXX(MXC), CYY(MXC)
      LOGICAL INVERS

      JJ = 1
 2    CONTINUE

        XMI = 1.D10
        DO 1 I = 1, NC
          IF(CX(I) .LT. XMI) THEN
            IMI = I
            XMI = CX(I)
          ENDIF
 1      CONTINUE

        CXX(JJ) = CX(IMI)
        CYY(JJ) = CY(IMI)

        IF(JJ .EQ. NC) GOTO 10
        CX(IMI) = 1.D10
        JJ = JJ+1
       
      GOTO 2

 10   CONTINUE

      INVERS = CYY(NC) .LT. CYY(1)

      XMI = 1.D12
      XMA = -1.D12
      YMI = 1.D12
      YMA = -1.D12
      IF(INVERS)  THEN
        II = NC
      ELSE
        II = 1
      ENDIF
      DO 3 I=1, NC
        IF(INVERS) THEN
          CX(I) =  - CXX(II)
COSY          CX(I) =  - CXX(II) / 10.D0
          CY(I) = CYY(II)
          II = II - 1
        ELSE
          CX(I) =  CXX(II)
COSY          CX(I) =  CXX(II) / 10.D0
          CY(I) = CYY(II)
          II = II + 1
        ENDIF

        IF(CX(I) .LT. XMI) XMI = CX(I)
        IF(CX(I) .GT. XMA) XMA = CX(I)
        IF(CY(I) .LT. YMI) YMI = CY(I)
        IF(CY(I) .GT. YMA) YMA = CY(I)
CCCCC        WRITE(*,*) I, CX(I), CY(I)
 3    CONTINUE

C      IF(ABS(YMI) .GT. ABS(YMA) ) THEN
C        SIGN =  -1.D0
C      ELSE
C        SIGN = 1.D0
C      ENDIF
          
C----- Normalise CY to 1
C and calculate approximate position of mag. face for fringe field application
      YN = YMA
C      IF(SIGN .EQ. -1.D0) THEN
C        YMA = -YMI
C        YMI = -YN
C        YN = YMA
C      ENDIF
        
      IF(NRM .EQ. 0) YN = 1.D0
      SUM = 0.D0        
      DX = CX(2) - CX(1)
      DO 4 I = 1, NC
        CY(I) = CY(I) / YN
        SUM = SUM + CY(I) * DX
 4    CONTINUE
      SUM = SUM - ( CY(1)+CY(NC) ) * DX/2.D0
      XFB = CX(NC) - SUM

      RETURN
      END
