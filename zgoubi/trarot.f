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
      SUBROUTINE TRAROT(TX,TY,TZ,RX,RY,RZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -------------------------------------------------
C     CHANGEMENT DE REFERENCE DE L'ENSEMBLE DU FAISCEAU
C     -------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YL2,ZL2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXSTEP.H'
      INCLUDE 'CSR.H'
C      COMMON/CSR/ KTRA,KCSR,YZXB(MXSTEP,41,36),DWC(MXT)
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
C     >,AMS  ,AMP,ENSTAR,BSTAR,TDVM,TETPHI(2,MXT)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      COMMON/DON/ A(09876,99),IQ(09876),IP(09876),NB,NOEL
      INCLUDE "MAXCOO.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      COMMON/GASC/ AI, DEN, KGA
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
 
      LOGICAL EVNT
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,100) TX,TY,TZ,RX,RY,RZ
 100    FORMAT(/,10X,' CHANGEMENT  DE  REFERENCIEL'
     >          ,1P,/,15X,' TRANSLATION:'
     >             ,/,15X,'   TX =',G12.4,'  m'
     >             ,/,15X,'   TY =',G12.4,'  m'
     >             ,/,15X,'   TZ =',G12.4,'  m'
     >             ,/,15X,' ROTATION:'
     >             ,/,15X,'   RX =',G12.4,'  rad'
     >             ,/,15X,'   RY =',G12.4,'  rad'
     >             ,/,15X,'   RZ =',G12.4,'  rad')
      ENDIF

      EVNT = KSPN .EQ. 1 .OR. IFDES .EQ. 1 .OR. KGA .EQ. 1 .OR.
     >  LIMIT .EQ. 1 .OR. KCSR.EQ. 1
 
      DO 1 IT=1,IMAX
 
C------- IEX<-1<=> PARTICULE STOPPEE
        IF( IEX(IT) .LT. -1) GOTO 1
 
        IF(IT .EQ. IREP(IT) .OR. .NOT.ZSYM) THEN
          CALL INITRA(IT)
          CALL TRROTE(EVNT,TX*100.D0,TY*100.D0,TZ*100.D0,RX,RY,RZ)
          CALL MAJTRA(IT)
        ELSE
          CALL DEJACA(IT)
        ENDIF
 
   1  CONTINUE
 
      IF(NRES.GT.0) WRITE(NRES,101) IEX(1),(F(J,1),J=1,7)
  101 FORMAT(' TRAJ 1 IEX,D,Y,T,Z,P,S,time :',I3,1P,5G12.4,2G17.5)

      RETURN
      END
