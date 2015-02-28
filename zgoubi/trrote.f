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
C  -------
      SUBROUTINE TRROTE(EVNT,TX,TY,TZ,RX,RY,RZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EVNT
C     --------------------------------------------------
C     TRANSLATION, ROTATION + TRANSPORT TO THE NEW FRAME
C     --------------------------------------------------
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CONST_3.H"      ! COMMON/CONST/ CL9,CEL,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
C      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
C     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
C     >,AMS ,AMP,ENSTAR,BSTAR,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
C      COMMON/GASC/ AI, DEN, KGA
      COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
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
c        if(kspn .eq. 1) then
c          sf(2,it) =   sf(2,it) * co + sf(3,it) * sn
c          sf(3,it) = - sf(2,it) * sn + sf(3,it) * co
c        endif
        IF(KSPN .EQ. 1 ) CALL SPNROT(IT,rx,ZERO,ZERO)
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
        if(kspn .eq. 1) 
     >  call endjob('SBR TRROTE. Spin rotation to be implemented. ',-99)
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
        if(kspn .eq. 1) 
     >  call endjob('SBR TRROTE. Spin rotation to be implemented. ',-99)
      ENDIF

C----- TRANSPORT TO NEW FRAME (=> X=0)
       Y=Y - X*TAN(T)
       Z=Z - X*TAN(P)/COS(T)
       SAR= SAR - X/(COS(T)*COS(P))

      IF(EVNT) CALL EVENT(-X,Y,T,Z,P,ZERO,UN,QBR,SAR,TAR,KEX,IT,
     > AMT,Q,BORO,I2,IMAX,*99)
C     > AMT,Q,BORO,I2,IFDES,KGA,I0,IMAX,*99)

 99   RETURN
      END
