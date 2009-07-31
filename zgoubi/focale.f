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
      SUBROUTINE FOCALE(IENERG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ***********************************************
C     RECHERCHE LES DIMENSIONS DU FAISCEAU EN XI AVEC
C     - POSITION XI DU WAIST CALCULEE
C          ET WAIST GLOBAL           SI IENERG=1
C          OU WAIST PAR MOMENTUM     SI IENERG=2
C     - LECTURE DE XI ET DIM. BEAM   SI IENERG=3
C     ***********************************************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT),AMQLU
      COMMON/FOCAL/ TTI(MXT),YI(MXT),ZI(MXT),WC,XI,YIO
     >,YMI,WCZ,MZ,IMAX1,IMAX2,MY
      COMMON /INIT/ FA0(6,6),FA1(6,6),BID(6),BID1(6),IF
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
 
      DIMENSION F2(MXT),F3(MXT),F4(MXT),F5(MXT)
      CHARACTER HOVE(2)*10, YZ(2)*1, LAHA(2)*7
      INTEGER HV
 
      DATA (HOVE(I),I=1,2) / 'HORIZONTAL' , 'VERTICAL' /
      DATA (YZ(I),I=1,2) / 'Y', 'Z' /
      DATA (LAHA(I),I=1,2) / 'LARGEUR' , 'HAUTEUR' /

      IF    (IENERG .GT.0) THEN
        HV=1
        DO 4 I=1,IMAX
          IF( IEX(I) .LT. -1 ) GOTO 4
          F2(I)=F(2,I)
          F3(I)=F(3,I)*.001D0
          F4(I)=F(4,I)
          F5(I)=F(5,I)*.001D0
 4      CONTINUE
        IEN=IENERG
      ELSEIF(IENERG .LT.0) THEN
        HV=2
        DO 5 I=1,IMAX
          IF( IEX(I) .LT. -1) GOTO 5
          F2(I)=F(4,I)
          F3(I)=F(5,I)*.001D0
          F4(I)=F(2,I)
          F5(I)=F(3,I)*.001D0
 5      CONTINUE
        IEN=-IENERG
      ENDIF
 
      IF    (IEN .EQ. 1 .OR. IEN .EQ. 3) THEN
        JDMAX=1
        JMAXT=IMAX
      ELSEIF(IEN .EQ. 2) THEN
        JDMAX=IDMAX
        JMAXT=IMAXT
      ENDIF
 
      DO 3 ID=1,JDMAX
        IMAX1=1+(ID-1)*JMAXT
        IMAX2=IMAX1+JMAXT-1
        STY=0D0
        ST2=0D0
        ST=0D0
        SY=0D0
        IMAXI=0
        DO 1 I=IMAX1,IMAX2
          IF( IEX(I) .LT. -1) GOTO 1
          TTI(I)=TAN(F3(I))
          STY=STY+TTI(I)*F2(I)
          ST2=ST2+TTI(I)*TTI(I)
          ST=ST+TTI(I)
          SY=SY+F2(I)
          IMAXI=IMAXI+1
 1      CONTINUE
 
        IF    (IEN .LT. 3) THEN
          XI=-(STY-(ST*SY)/IMAXI)/(ST2-(ST*ST)/IMAXI)
        ELSEIF(IEN .EQ. 3) THEN
          XI = A(NOEL,1)
        ENDIF
 
        YMI=0D0
        YMIN= 1.D10
        YMAX=-1.D10
        ZMAX=0D0
        WI=0D0
        YIO=F2(IMAX1)+XI*TTI(IMAX1)
        DO 2 I=IMAX1,IMAX2
          IF( IEX(I) .LT. -1 ) GOTO 2
          YI(I)=F2(I)+XI*TTI(I)-YIO
          ZI(I)=F4(I)+XI*TAN(F5(I))
C          SARI(I)= F(6,I)+XI/( COS(F3(I))*COS(F5(I)) )
          YMI=YMI+YI(I)
          IF(YI(I) .LT. YMIN) YMIN=YI(I)
          IF(YI(I) .GT. YMAX) YMAX=YI(I)
C          ZMAX=AMAX1(ZMAX,ABS(ZI(I)))
          IF(ZMAX .LT. ABS(ZI(I))) ZMAX=ABS(ZI(I))
          WI=WI+YI(I)*YI(I)
 2      CONTINUE
        YMI=YMI/(DBLE(IMAXI))
        WI=   SQRT((WI/DBLE(IMAXI))-YMI*YMI)*2.35D0
        WC=0.005D0
C       IF(WI.GT.0.1D0) WC=0.01*DBLE(IFIX(WI*10.D0))*0.5
        IF(WI.GT.0.1D0) WC=0.01D0*INT(WI*10.D0)*0.5D0
        MY=200.D0*WC+0.5D0
C       WCZ=0.05*DBLE(IFIX(ZMAX/0.05)/30+1)
        WCZ=0.05D0*(NINT(ZMAX/0.05D0)/30+1)
        MZ=100.D0*WCZ+0.5D0
        WT=YMAX-YMIN
        SWI2 = SWI2 + WI*WI
 
        IF(NRES.GT.0) THEN
          WRITE(NRES,103)HOVE(HV),IMAXI,JMAXT,HOVE(HV)
     >    ,XI,YZ(HV),YIO,YZ(HV),YMI,LAHA(HV),WI,WT
          CALL IMPTRA(IMAX1,IMAX2,NRES)
C          CALL TRACE
        ENDIF
 
    3 CONTINUE
 
      RETURN
 
  103 FORMAT(/,1P,
     >5X,'RECHERCHE DU POINT DE FOCALISATION ',A,' DE ',I6
     1,' TRAJECTOIRES (SUR ',I6,')'
     1,//,5X,'POINT DE FOCALISATION ',A
     2,' SUR L ORBITE MOYENNE    X =',G14.6,' CM    ',A,' =',G14.6
     3,' CM',//,5X,'DECALAGE DU CENTRE DE GRAVITE EN ',A,' = ',G14.6
     4,' CM',//,5X,A,' IMAGE,  A MI-HAUTEUR =',G14.6,' CM,  TOTALE ='
     5,G14.6,' CM',/,21X,2(17X,9('.')),/)
 
      END
