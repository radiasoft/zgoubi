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
      SUBROUTINE AUTORF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------------------------
C     CHANGEMENT AUTOMTIQ DE REFERENTIEL
C     ------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT),AMQLU
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
 
      IOP = A(NOEL,1)
      IF(IOP .EQ. 3) THEN
        IRF=A(NOEL,10)
        MX1=A(NOEL,11)
        MX2=A(NOEL,12)
      ELSE
        IRF=1
      ENDIF
 
      IF    (IOP .EQ. 1) THEN
C------- POSITIONNEMENT REFERENCE = TRAJ. 1, ABSCISSE ACTUELLE
        XC =ZERO
        YC =F(2,IRF)
      ELSEIF(IOP .EQ. 2) THEN
C------- POSITIONNEMENT REFERENCE = WAIST DES TRAJ. 1 ET 4-5
        CALL FOCAL1(IRF,4,5,XC,YC )
      ELSEIF(IOP .EQ. 3) THEN
C------- POSITIONNEMENT REFERENCE = WAIST DES TRAJ. IRF ET MX1-MX2
        CALL FOCAL1(IRF,MX1,MX2,XC,YC )
      ENDIF
      AA  =F(3,IRF)*.001D0
 
      DO 1 I=1,IMAX
C       +++ IEX<-1 <=> Particule stoppee
        IF( IEX(I) .LT. -1) GOTO 1
        IF(I .EQ. IREP(I) .OR. .NOT.ZSYM) THEN
          CALL INITRA(I)
          CALL CHAREF(.FALSE.,XC,YC,AA)
          CALL MAJTRA(I)
        ELSE
          CALL DEJACA(I)
        ENDIF
   1  CONTINUE
 
      IF(NRES .GT. 0) THEN
        WRITE(NRES,100) XC,YC,AA*DEG,AA
 100    FORMAT(/,' Change  of  reference   XC =',F9.3,' cm , YC =',
     >   F10.3,' cm ,   A =',F12.5,' deg  (i.e., ',F10.6,' rad)',/)
C 100    FORMAT(/,' CHANGEMENT  DE  REFERENCE  XC =',F9.3,' cm , YC =',
C     >   F10.3,' cm ,   A =',F12.5,' deg  (i.e., ',F10.6,' rad)',/)
        WRITE(NRES,101) IEX(1),(F(J,1),J=1,7)
  101   FORMAT(' TRAJ 1 IEX,D,Y,T,Z,P,S,time :',I3,1P,5G12.4,2G17.5)
      ENDIF
 
      RETURN
      END
