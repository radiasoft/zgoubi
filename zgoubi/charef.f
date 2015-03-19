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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE CHAREF(EVNT,XC,YC,A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EVNT
C     -----------------------------------------------
C     CHANGEMENT DE REFERENCE PARTICULE PAR PARTICULE
C     -----------------------------------------------
      INCLUDE "C.AIM_2.H"     ! COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE "MAXTRA.H"
      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
C     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      COMMON/GASC/ AI, DEN, KGA
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      INCLUDE "C.TRAJ.H"     ! COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      PARAMETER(I0=0)
      YO=Y
      Y=((Y-YC)*COS(T)+XC*SIN(T))/COS(T-A)
      T=T-A 
      XL=XC-Y*SIN(A)
      YL=YC-YO+Y*COS(A)
      DL=SQRT(XL*XL+YL*YL)
      DL=SIGN(DL,XL)
      DS = DL/COS(P)
      SAR= SAR+DS
      Z=Z+DL*TAN(P)
      QBRO = QBR*CL9
      DTAR = DS / (QBRO/SQRT(QBRO*QBRO+AMT*AMT)*CL9)
      TAR = TAR + DTAR
      IF(EVNT) THEN
C------- spin tracking
        IF(KSPN .EQ. 1 ) THEN
          IF(A .NE. 0.D0) CALL SPNROT(IT,ZERO,ZERO,A)
        ENDIF
        YY = Y
        XX = X
        IF(KART.EQ.2) THEN
C--------- Cylindrical coordinates
          YY = Y+RM
          XX = A
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
