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
      SUBROUTINE CHANRF(NSR,EVNT,QSHRO,VSHRO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EVNT
      CHARACTER QSHRO(*)*(*)
      DIMENSION VSHRO(*)
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

C        WRITE(*,FMT='(A)') ' SBR CHANRF ',(qshro(I),I=1,NSR),NSR
C              STOP

      ksr = 1
 1    continue
        if    (qshro(ksr).eq.'XS') then
c          write(77,*) 'chanrf noel,ksr,qshro(ksr): ',
c     >       noel,ksr,qshro(ksr),VSHRO(ksr),nsr
          XS = VSHRO(ksr) 
          YS=0.D0
          ZR = 0.D0

          YO=Y
          Y=((Y-YS)*COS(T)+XS*SIN(T))/COS(T-ZR)
          T=T-ZR 
          XL=XS-Y*SIN(ZR)
          YL=YS-YO+Y*COS(ZR)
          DL=SQRT(XL*XL+YL*YL)
          DL=SIGN(DL,XL)
          DS = DL/COS(P)
          SAR= SAR+DS
          Z=Z+DL*TAN(P)
          QBRO = QBR*CL9
          DTAR = DS / (QBRO/SQRT(QBRO*QBRO+AMT*AMT)*CL9)
          TAR = TAR + DTAR
c          YO=Y
c          XL=XC
c          DL=SQRT(XL*XL)
c          DL=SIGN(DL,XL)
c          DS = DL/COS(P)
c          SAR= SAR+DS
c          Z=Z+DL*TAN(P)
c          QBRO = QBR*CL9
c          DTAR = DS / (QBRO/SQRT(QBRO*QBRO+AMT*AMT)*CL9)
c          TAR = TAR + DTAR
        elseif(qshro(ksr).eq.'YS') then
c          write(77,*) 'chanrf noel,ksr,qshro(ksr): ',
c    >       noel,ksr,qshro(ksr),VSHRO(ksr),nsr
          XS = 0.D0
          YS = VSHRO(ksr) 
          ZR = 0.D0
          Y = Y-YS
        elseif(qshro(ksr).eq.'ZS') then
          ZC = VSHRO(ksr) 
          Z = Z-ZC
        elseif(qshro(ksr).eq.'XR') then
          XR = VSHRO(ksr) 
        elseif(qshro(ksr).eq.'YR') then
          YR = VSHRO(ksr) 
        elseif(qshro(ksr).eq.'ZR') then
          write(77,*) 'chanrf noel,ksr,qshro(ksr): ',
     >       noel,ksr,qshro(ksr),VSHRO(ksr),nsr 
          XS = 0.D0
          YS = 0.D0
          ZR = VSHRO(ksr) 

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
        endif

        KSR = KSR+1
        IF(KSR.LE.NSR) GOTO 1

      IF(EVNT) THEN
C------- spin tracking
        IF(KSPN .EQ. 1 ) THEN
          IF(ZR .NE. 0.D0) CALL SPNROT(IT,ZERO,ZERO,ZR)
        ENDIF
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
      ENDIF
C 99   RETURN 
      RETURN 
      END
