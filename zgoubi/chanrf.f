C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory   
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
C      SUBROUTINE CHANRF(NSR,EVNT,QSHRO,VSHRO)
      SUBROUTINE CHANRF(EVNT,QSHRO,VSHRO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EVNT
      PARAMETER (MSR=10)
      CHARACTER QSHRO(MSR)*(2)
      DIMENSION VSHRO(MSR)
C     -----------------------------------------------
C     CHANGEMENT DE REFERENCE PARTICULE PAR PARTICULE
C     -----------------------------------------------
      COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE "MAXTRA.H"
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/GASC/ AI, DEN, KGA
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      PARAMETER(I0=0)

      NSR = NINT(VSHRO(MSR))

      KSR = 1
 1    CONTINUE
        IF    (QSHRO(KSR).EQ.'XS') THEN
          XS = VSHRO(KSR) 
          YO=Y
          Y=(Y*COS(T)+XS*SIN(T))/COS(T)
          XL=XS
          YL=-YO+Y
          DL=SQRT(XL*XL+YL*YL)
          DL=SIGN(DL,XL)
          DS = DL/COS(P)
          SAR= SAR+DS
          Z=Z+DL*TAN(P)
          QBRO = QBR*CL9
          DTAR = DS / (QBRO/SQRT(QBRO*QBRO+AMT*AMT)*CL9)
          TAR = TAR + DTAR
        ELSEIF(QSHRO(KSR).EQ.'YS') THEN
          YS = VSHRO(KSR) 
          Y = Y-YS
        ELSEIF(QSHRO(KSR).EQ.'ZS') THEN
          ZS = VSHRO(KSR) 
          Z = Z-ZS
        ELSEIF(QSHRO(KSR).EQ.'XR') THEN
          XR = VSHRO(KSR) 
          CALL ENDJOB('SBR CHANRF. XR is not implemented.',-99)
C          IF(KSPN .EQ. 1 ) CALL SPNROT(IT,+-XR,ZERO,ZERO)
        ELSEIF(QSHRO(KSR).EQ.'YR') THEN
          YR = VSHRO(KSR) 
          ZO=Z
C          pp = atan2(tan(P),cos(T))
          pp = atan(tan(P)/cos(T))
          Z=Z*COS(pp)/COS(pp-YR)
C          pp=pp-YR 
C          P =  atan2(tan(pp),1.d0/cos(T))
          XL = -Z *SIN(YR)
          ZL = -ZO + Z*COS(YR)
          DL=SQRT(XL*XL+ZL*ZL)
          DL=SIGN(DL,XL)
C          DS = DL/sin(pp)/SIN(P)
          DS = DL/COS(T)/COS(P)
          SAR= SAR+DS
          Y=Y+DS*COS(P)*SIN(T)
          QBRO = QBR*CL9
          DTAR = DS / (QBRO/SQRT(QBRO*QBRO+AMT*AMT)*CL9)
          TAR = TAR + DTAR
C          CALL ENDJOB('SBR CHANRF. YR is not implemented.',-99)
C          IF(KSPN .EQ. 1 ) CALL SPNROT(IT,ZERO,+-YR,ZERO)
        ELSEIF(QSHRO(KSR).EQ.'ZR') THEN
          ZR = VSHRO(KSR) 
          YO=Y
          Y=Y*COS(T)/COS(T-ZR)
          T=T-ZR 
          XL=-Y*SIN(ZR)
          YL=-YO+Y*COS(ZR)
          DL=SQRT(XL*XL+YL*YL)
          DL=SIGN(DL,XL)
          DS = DL/COS(P)
          SAR= SAR+DS
          Z=Z+DL*TAN(P)
          QBRO = QBR*CL9
          DTAR = DS / (QBRO/SQRT(QBRO*QBRO+AMT*AMT)*CL9)
          TAR = TAR + DTAR
          IF(KSPN .EQ. 1 ) CALL SPNROT(IT,ZERO,ZERO,ZR)
        ENDIF

        KSR = KSR+1
      IF(KSR.LE.NSR) GOTO 1

      IF(EVNT) THEN

CC------- spin tracking
C        IF(KSPN .EQ. 1 ) THEN
C          IF(ZR .NE. 0.D0) CALL SPNROT(IT,ZERO,ZERO,ZR)
C        ENDIF
        YY = Y
        XX = X
        IF(KART.EQ.2) THEN
C--------- Cylindrical coordinates
          YY = Y + RM
          XX = ZR
        ENDIF

C Problems with DIPOLE-M when calling EVENT/CHAMBRE here : y itself can be 
C either y or y+rm  depending when it is called
C            write(*,*) ' sbr charef ',it,y,rm,y+rm
C        CALL EVENT(DL,YY,T,Z,P,XX,UN,QBR,SAR,TAR,KEX,IT,
C     >  AMT,QT,BORO,KART,IFDES,KGA,I0,IMAX,*99)
C 99   CONTINUE

      ENDIF

      RETURN 
      END
