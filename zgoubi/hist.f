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
C  Brookhaven National Laboratory               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE HIST
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------------------
C     CONSTITUE UN TABLEAU POUR LISTING
C     LE TABLEAU EST LISTE PAR HISTO
C     ---------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,40)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      PARAMETER (JH=24,KH=5)
      COMMON/HISTO/ ICTOT(JH,KH),MOYC(JH,KH) ,CMOY(JH,KH),JMAX(JH,KH)
      COMMON/HISTOG/ NC(JH,120,KH),J,NH,XMI(JH,KH),XMO(JH,KH),XMA(JH,KH)
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IMAXD,IMAXT
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
 
      DIMENSION XMO2(JH,KH), IMX(JH,KH)
      CHARACTER TYP(JH,KH)*2 ,TEXT*38 ,KAR
      CHARACTER NORME(2)*14
      CHARACTER KOORD(JH)*6
      CHARACTER KUNIT(JH)*6
 
      DATA NORME /'NORMALISE     ','NON  NORMALISE'/
      DATA KOORD /'  D  ','Y    ','THETA','Z    ','PHI  ','  S  '
     >,4*' Time ','  Do ','Yo   ','To   ','Zo   ','PHIo ','  So '
     >,4*' Time ','Sx   ','Sy   ','Sz   ','<S>  ' /
      DATA KUNIT /'      ',' (CM) ',' (MRD)',' (CM) ',' (MRD)',' (CM) '
     >,4*' (s)  ','      ',' (CM) ',' (MRD)',' (CM) ',' (MRD)',' (CM) '
     >,8*' (s)  ' /
      DATA MIDCOL/ 50 /
 
 
 
C                   NUM.COORD.  BORNES      NCOL/LISTING     NUMERO DE
C                     1-JH     PHYSIQUES      < 120           L'HISTO
C                              COMPTAGE                        1-KH
C      READ(NDAT,*)    J     , XMIN,XMAX,      NCOL        ,    NH
      J = A(NOEL,1)
      XMIN = A(NOEL,2)
      XMAX = A(NOEL,3)
      NCOL = A(NOEL,4)
      NH = A(NOEL,5)
C                                                     COMPTE LES S(SEC),
C                   AMPLI.VERT.  SYMBOLE  NORM.VERT   P(PRIM) OU Q (PAS
C                    #30                    1-2       DE SELECTION)
C      READ(NDAT,*) NBLINE    ,   KAR   ,  NORMY  ,       TYP(J,NH)
      NBLINE = A(NOEL,10)
      KAR = TA(NOEL,1)(1:1)
      NORMY = A(NOEL,11)
      TYP(J,NH) = TA(NOEL,2)(1:2)
C
C     LE GRAPHIQ SE CENTRE AUTOMATIQT SUR LA COLONNE
C     MIDCOL DU LSTING. SA LARGEUR EST DE NCOL COLONNES
C     POUR DES BORNES PHYSIQUES DE VALEURS XMIN ET XMAX
 
      IF(IPASS .EQ. 1) THEN
        IMX(J,NH)=0
        XMO2(J,NH)=0D0
      ENDIF
 
      IF(NCOL.GT.2*MIDCOL) NCOL=2*MIDCOL
      FNORM=NCOL/(XMAX-XMIN)
      NC2=NCOL/2
      IC1=MIDCOL-NC2+1
      IC2=IC1+NCOL
      JMAX(J,NH) = JMAX(J,NH) + IMAX
 
C     ** LIMITES ET MOYENNE PHYS. DE LA VARIABLE :
      DO 3 I=1,IMAX
 
C       +++ IEX<-1<=> PARTICULE STOPPEE
        IF(IEX(I) .LT. -1) GOTO 3
 
        IF    (J .LE.  6) THEN
C         ** DP/P initial, Y, T, Z, P, S
          FJI = F(J,I)
        ELSEIF(J .LE. 16) THEN
C         ** DP/P final Yo, To, Zo, Po, So
          FJI = FO(J-10,I)
        ELSEIF(J .LE. 24) THEN
C         ** COMPOSANTES SPIN: SX, SY, SZ, ET <S>=SQRT(SX2+SY2+SZ2)
          FJI = SF(J-20,I)
        ENDIF
 
C        IF    (IFDES .EQ. 1) THEN
          IF(  (TYP(J,NH) .EQ. 'P' .AND. LET(I) .NE. 'S')
     >    .OR. (TYP(J,NH) .EQ. 'S' .AND. LET(I) .EQ. 'S')
     >    .OR. (TYP(J,NH) .EQ. 'Q')  ) THEN
             IF(FJI .LT. XMI(J,NH)) XMI(J,NH) = FJI
             IF(FJI .GT. XMA(J,NH)) XMA(J,NH) = FJI
             IMX(J,NH) = IMX(J,NH) + 1
             XMO(J,NH) = XMO(J,NH) + FJI
             XMO2(J,NH) = XMO2(J,NH) + FJI*FJI
          ENDIF
C        ELSEIF(IFDES .EQ. 0) THEN
C          IF(FJI .LT. XMI(J,NH)) XMI(J,NH) = FJI
C          IF(FJI .GT. XMA(J,NH)) XMA(J,NH) = FJI
C          IMX(J,NH) = IMX(J,NH) + 1
C          XMO(J,NH) = XMO(J,NH) + FJI
C          XMO2(J,NH) = XMO2(J,NH) + FJI*FJI
C        ENDIF
 
 3    CONTINUE
 
C     ** COMPTAGE DANS LA FENETRE  XMIN-XMAX
      DO 4 I=1,IMAX
 
C       +++ IEX<-1<=> PARTICULE STOPPEE
        IF(IEX(I) .LT. -1) GOTO 4
 
        IF    (J .LE.  6) THEN
C         ** DP/P initial, Y, T, Z, P, S
          FJI = F(J,I)
        ELSEIF(J .LE. 16) THEN
C         ** DP/P final Yo, To, Zo, Po, So
          FJI = FO(J-10,I)
        ELSEIF(J .LE. 24) THEN
C         ** SPIN: SX, SY, SZ, ET <S>=SQRT(SX2+SY2+SZ2)
          FJI = SF(J-20,I)
        ENDIF
 
C        IF    (IFDES .EQ. 1) THEN
          IF(  (TYP(J,NH) .EQ. 'P' .AND. LET(I) .NE. 'S')
     >    .OR. (TYP(J,NH) .EQ. 'S' .AND. LET(I) .EQ. 'S')
     >    .OR. (TYP(J,NH) .EQ. 'Q')  ) THEN
             IF(FJI .GE. XMIN .AND. FJI .LE. XMAX) THEN
               ICOL = IC1 + NINT( ( FJI - XMIN ) * FNORM )
               NC(J,ICOL,NH) = NC(J,ICOL,NH) + 1
             ENDIF
          ENDIF
C        ELSEIF(IFDES .EQ. 0) THEN
C          IF(FJI .GE. XMIN .AND. FJI .LE. XMAX) THEN
C            ICOL = IC1 + NINT( ( FJI - XMIN ) * FNORM )
C            NC(J,ICOL,NH) = NC(J,ICOL,NH) + 1
C          ENDIF
C        ENDIF
 
 4    CONTINUE
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,106) KOORD(J)
 106    FORMAT(/,30X,'HISTOGRAMME  DE  LA  COORDONNEE  ',A)
C        IF    (IFDES .EQ. 1) THEN
          IF     (TYP(J,NH) .EQ. 'P') THEN
           TEXT='PARTICULES  PRIMAIRES'
          ELSEIF (TYP(J,NH) .EQ. 'S') THEN
           TEXT='PARTICULES  SECONDAIRES'
          ELSEIF (TYP(J,NH) .EQ. 'Q') THEN
           TEXT='PARTICULES  PRIMAIRES  ET  SECONDAIRES'
          ENDIF
          WRITE(NRES,104) TEXT
 104      FORMAT(30X,A)
C        ENDIF
 
        WRITE(NRES,103) XMIN,XMAX,KUNIT(J),NORME(NORMY)
 103    FORMAT(30X,'DANS  LA  FENETRE : '
     >  ,1P,G12.4,' / ',G12.4,A,/,30X,A,/)
 
        CALL HISTAB(IC1,IC2, NBLINE, KAR, NORMY, *5)
 
        POSIT=(MOYC(J,NH)-IC1)/FNORM + XMIN
        WRITE(NRES,102) POSIT,KUNIT(J), 1.D0/FNORM,KUNIT(J)
 102    FORMAT(15X,' VAL. PHYS. AU  "      "         : ',1P,G10.3,A6,/
     >        ,15X,' RESOLUTION  PAR  CANAL          : ',   G10.3,A6)
 
 5      CONTINUE
        SIGMA=SQRT(XMO2(J,NH)/IMX(J,NH)-(XMO(J,NH)/IMX(J,NH))**2)
        WRITE(NRES,105)
     >  IMX(J,NH),XMI(J,NH),XMA(J,NH),XMA(J,NH)-XMI(J,NH)
     >  ,KUNIT(J),XMO(J,NH)/IMX(J,NH), KUNIT(J), SIGMA,KUNIT(J)
 105    FORMAT(/,15X,' PARAMETRES  PHYSIQUES  DE  LA  DISTRIBUTION :'
     >  ,/,30X,   'COMPTAGE =',I7,'  PARTICULES'
     >  ,/,30X,1P,'MIN =',G12.4,', MAX =',G12.4,', MAX-MIN =',G12.4,A
     >  ,/,30X,   'MOYENNE =',G12.4,A
     >  ,/,30X,   'SIGMA =',G12.4,A,/)
        WRITE(NRES,101) IEX(1),(F(J,1),J=1,7)
  101   FORMAT(' TRAJ 1 IEX,D,Y,T,Z,P,S,time :',I3,1P,5G12.4,2G17.5)

      ENDIF
 
      RETURN
      END
