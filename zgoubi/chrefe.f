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
      SUBROUTINE CHREFE(IOP,NRES,XC,YC,AA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -------------------------------------------------
C     CHANGEMENT DE REFERENCE DE L'ENSEMBLE DU FAISCEAU
C     -------------------------------------------------
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YL2,ZL2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXSTEP.H'
      INCLUDE 'CSR.H'
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      COMMON/GASC/ AI, DEN, KGA
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      COMMON/SYNRA/ KSYN
 
      LOGICAL EVNT
 
      EVNT = KSPN.EQ.1 .OR. IFDES.EQ.1 .OR. KGA.EQ.1 .OR. 
     >  LIMIT.EQ.1 .OR. KSYN.GE.1 .OR. KCSR.EQ.1 

      DO 1 IT=1,IMAX
C------- IEX<-1<=> PARTICULE STOPPEE
        IF( IEX(IT) .LT. -1) GOTO 1
 
        IF(IT .EQ. IREP(IT) .OR. .NOT.ZSYM) THEN
          CALL INITRA(IT)
          CALL CHAREF(EVNT,XC,YC,AA)
          CALL MAJTRA(IT)
        ELSE
          CALL DEJACA(IT)
        ENDIF
 
   1  CONTINUE
 
      IF(IOP .EQ. 1) THEN
        IF(NRES.GT.0) THEN
          WRITE(NRES,100) XC,YC,AA*DEG,AA
 100      FORMAT(/,' CHANGEMENT  DE  REFERENCE  XC ='
     >    ,F10.3,' CM , YC ='
     >    ,F10.3,'  CM ,   A =',F12.5,' DEG  (i.e.,',F10.6,' rad)',/)
          WRITE(NRES,101) IEX(1),(F(J,1),J=1,7)
  101     FORMAT(' TRAJ 1 IEX,D,Y,T,Z,P,S,time :',I3,1P,5G12.4,2G17.5)
        ENDIF
      ENDIF
 
      RETURN
      END
