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
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
 
      IOP = NINT(A(NOEL,1))
C      IF    (IOP .EQ. 3) THEN
C        IRF=NINT(A(NOEL,10))
C        MX1=NINT(A(NOEL,11))
C        MX2=NINT(A(NOEL,12))
C      ELSE
C        IRF=1
C      ENDIF
 
      IF    (IOP .LE. 3) THEN
        IF    (IOP .EQ. 0) THEN
          GOTO 99
        ELSEIF(IOP .EQ. 1) THEN
C--------- POSITIONNEMENT REFERENCE = TRAJ. 1, ABSCISSE ACTUELLE
          IRF=1
          XC =ZERO
          YC =F(2,IRF) 
        ELSEIF(IOP .EQ. 2) THEN
C--------- POSITIONNEMENT REFERENCE = WAIST DES TRAJ. 1 ET 4-5
          IRF=1
          CALL FOCAL1(IRF,4,5,
     >                        XC,YC )
        ELSEIF(IOP .EQ. 3) THEN
C--------- POSITIONNEMENT REFERENCE = WAIST DES TRAJ. IRF ET MX1-MX2
          IRF=NINT(A(NOEL,10))
          MX1=NINT(A(NOEL,11))
          MX2=NINT(A(NOEL,12))
          CALL FOCAL1(IRF,MX1,MX2,
     >                            XC,YC )
        ENDIF
        AA  =F(3,IRF) * 0.001D0
      ELSEIF(IOP .EQ. 4) THEN
        XC = ZERO
        YC = ZERO
        AA = ZERO
        II = 0
        DO I = 1, IMAX
          IF( IEX(I) .GT. 0) THEN
            II = II + 1
            YC = YC + F(2,I) 
            AA = AA + F(3,I) 
          ENDIF
        ENDDO
        YC = YC / DBLE(II)        
        AA = AA / DBLE(II) * 0.001D0 
      ELSE
        CALL ENDJOB('Sbr autorf. No such option I = ',IOP)
      ENDIF

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
 
 99   CONTINUE

      IF(NRES .GT. 0) THEN

        IF(IOP .EQ. 0) THEN
          WRITE(NRES,102) 
 102      FORMAT(/,20X,' AUTOREF  is  off')
        ELSE
          WRITE(NRES,100) XC,YC,AA*DEG,AA
 100      FORMAT(/,' Change  of  reference   XC =',F15.8,' cm , YC =',
     >    F15.8,' cm ,   A =',F16.9,' deg  (i.e., ',F14.10,' rad)',/)
C 100      FORMAT(/,' CHANGEMENT  DE  REFERENCE  XC =',F9.3,' cm , YC =',
C     >     F10.3,' cm ,   A =',F12.5,' deg  (i.e., ',F10.6,' rad)',/)
          WRITE(NRES,101) IEX(1),(F(J,1),J=1,7)
  101     FORMAT(' TRAJ 1 IEX,D,Y,T,Z,P,S,time :',I3,1P,5G12.4,2G17.5)
        ENDIF
      ENDIF
 
      RETURN
      END
