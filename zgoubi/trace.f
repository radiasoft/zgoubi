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
      SUBROUTINE TRACE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT),AMQLU
      COMMON/FOCAL/TTI(MXT),YI(MXT),ZI(MXT),WC,XI,YIO,YMI
     A,WCZ,MZ,IMAX1,IMAX2,MY
C      COMMON/LIGNE/ NT(200),NC(200)
      COMMON/LIGNE/ NT(MXT),NC(MXT)
 
      WRITE(NRES,100) MY,MZ,XI,YIO,YMI
 
      YM=65.5D0*WC
      DO 1 L=1,61
        ZLP=(31.5D0-L)*WCZ
        ZLM=ZLP-WCZ
        NKAR=0
        DO 2 I=IMAX1,IMAX2
          IF( IEX(I) .LT. -1) GOTO 2
          IF(ZI(I).GT.ZLP.OR.ZI(I).LE.ZLM) GOTO 2
          NK=((YM-YI(I))/WC+1.D0)
          IF(NK.LT.1.OR.NK.GT.131)  GOTO 2
          NKAR=NKAR+1
          IF(NKAR.GT.MXT) 
     >      STOP " *** Error, SBR TRACE -> NKAR should be < MXT"
          NT(NKAR)=I
          NC(NKAR)=NK
    2   CONTINUE
        K=2
        IF(MOD(L,10) .EQ. 1)  K=1
        CALL LIGN(K,NKAR)
    1 CONTINUE
 
      WRITE(NRES,101) MY,XI,YIO,YMI
 
      DO 3 L=1,61
        NKAR=0
        K=2
        IF(MOD(L,6) .NE. 1)  GOTO 4
        IF(MOD(L,30) .EQ. 1)  K=1
        DX=(31-L)*.3333333333D0
        DO 5 I=IMAX1,IMAX2
          IF( IEX(I) .LT. -1) GOTO 5
          IF(F(5,I) .NE. 0.D0)  GOTO 5
          YT=YI(I)+DX*TTI(I)
          NK=((YM-YT   )/WC+1.D0)
          IF(NK.LT.1.OR.NK.GT.131)  GOTO 5
          NKAR=NKAR+1
          NT(NKAR)=I
          NC(NKAR)=NK
    5   CONTINUE
    4 CONTINUE
      CALL LIGN(K,NKAR)
    3 CONTINUE
 
      RETURN
 
  100 FORMAT(1H1,3X,' ABERRATIONS VERTICALES',/
     1,3X,' MAILLE EN Y=' ,I4,' MM     EN Z ='
     1,I4,' MM   POINT CENTRAL X=',1P,E8.2,' CM Y=',E8.2,' CM',/
     2,3X,' CENTRE DE GRAVITE DECALE DE',E8.2,' CM')
  101 FORMAT(1H1,3X,' TRAJECTOIRES DU PLAN MEDIAN',/
     1,3X,' MAILLE EN Y=' ,I4
     1,' MM   EN X=10CM    POINT CENTRAL  X=',1PE8.2,' CM  Y=',E8.2
     1,' CM',/,3X,' CENTRE DE GRAVITE DECALE DE',E8.2,' CM')
      END
