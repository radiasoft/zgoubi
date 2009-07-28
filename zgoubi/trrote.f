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
      SUBROUTINE TRROTE(EVNT,TX,TY,TZ,RX,RY,RZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EVNT
C     --------------------------------------------------
C     TRANSLATION, ROTATION + TRANSPORT TO THE NEW FRAME
C     --------------------------------------------------
      INCLUDE "MAXTRA.H"
      COMMON/CONST/ CL9,CEL,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
C     >,AMS ,AMP,ENSTAR,BSTAR,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      COMMON/GASC/ AI, DEN, KGA
      COMMON/RIGID/ BORO,DPREF,DP,BR
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      PARAMETER(MT=41,MP=41,MO=101)
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      PARAMETER (I0=0,I2=2)

C----- TRANSLATION
      X=-TX
      Y=Y-TY
      Z=Z-TZ

C----- ROTATION
      IF(RX .NE. 0.D0) THEN
        CP=COS(P)
        UX=CP*COS(T)
        UY=CP*SIN(T)
        UZ=SIN(P)         
        SN=SIN(RX)
        CO=COS(RX)
        Y1=Y
        Y=  Y*CO+Z*SN
        Z=-Y1*SN+Z*CO
        U1=UY
        UY= UY*CO+UZ*SN
        UZ=-U1*SN+UZ*CO
        T=ATAN2(UY,UX)
        P=ATAN( UZ/SQRT(UX*UX+UY*UY) )
      ENDIF
      IF(RY .NE. 0.D0) THEN
        CP=COS(P)
        UX=CP*COS(T)
        UY=CP*SIN(T)
        UZ=SIN(P)         
        SN=SIN(RY)
        CO=COS(RY)
        X1=X
        X=  X*CO-Z*SN
        Z= X1*SN+Z*CO
        U1=UX
        UX= UX*CO-UZ*SN
        UZ= U1*SN+UZ*CO
        T=ATAN2(UY,UX)
        P=ATAN( UZ/SQRT(UX*UX+UY*UY) )
      ENDIF
      IF(RZ .NE. 0.D0) THEN
        CP=COS(P)
        UX=CP*COS(T)
        UY=CP*SIN(T)
        UZ=SIN(P)         
        SN=SIN(RZ)
        CO=COS(RZ)
        X1=X
        X=  X*CO+Y*SN
        Y=-X1*SN+Y*CO
        U1=UX
        UX= UX*CO+UY*SN
        UY=-U1*SN+UY*CO
        T=ATAN2(UY,UX)
        P=ATAN( UZ/SQRT(UX*UX+UY*UY) )
      ENDIF

C----- TRANSPORT TO NEW FRAME (=> X=0)
       Y=Y - X*TAN(T)
       Z=Z - X*TAN(P)/COS(T)
       SAR= SAR - X/(COS(T)*COS(P))

      IF(EVNT) CALL EVENT(-X,Y,T,Z,P,ZERO,UN,BR,SAR,TAR,KEX,IT,
     > AMT,QT,BORO,I2,IFDES,KGA,I0,IMAX,*99)

 99   RETURN
      END
