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
      SUBROUTINE DIPI(SCAL, 
     >                     DSREF)
C     >                      XL,DEV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C---------------------------------------------------------
C     Same as DIPSI, but for a single dipole. 
C---------------------------------------------------------
      COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
C      INCLUDE "MAXTRA.H"
C      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,HNORM
C     > ,YCH,ZCH
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DPPP,QBR,BRI
 
      PARAMETER (NMAG=5)
      DIMENSION CE(NMAG,6),CS(NMAG,6),C3(NMAG,6)
 
      DIMENSION FTAB(5,5)
      
      DOUBLE PRECISION LAMBDE, LAMBDS, LAMBD3
 
      SAVE CE,CS,C3,SHIFTE,SHIFTS,SHIFT3
      SAVE COEFN,COEFB,COEFG,ACN
      SAVE LAMBDE,QSIE,NCOEFE
      SAVE LAMBDS,QSIS,NCOEFS
      SAVE LAMBD3,QSI3,NCOEF3
      SAVE UMEGA,THETA,R1,U1,U2,R2
      SAVE UMEGAS,THETAS,R1S,U1S,U2S,R2S
      SAVE NN, RESOL

      LOGICAL OKTTA, LTXI, GTXF
      DATA OKTTA, LTXI, GTXF / .TRUE., .FALSE.   , .FALSE.   /

      NP = 2
      AT    = A(NOEL,NP)
      NP=NP+1 
      RM    = A(NOEL,NP)

C  HNORM=Champ MAX DANS LE DIPOLE.
C  COEFN=N=INDICE DE Champ, B=N', G=N''.
C  AT=ANGLE TOTAL DE LA CARTE DE Champ, ACENT='ANGLE AU CENTRE',
C    RM,MIN,MAX=RAYONS MOYEN,MIN,MAX DE LACARTE DE Champ. 
C  NBFACE=(2)3 : dipole limited by (2)3 field boundaries

      KMAG = 1

      NP=NP+1 
      ACENT = A(NOEL,NP)
      NP=NP+1 
      HNORM = A(NOEL,NP)*SCAL
      NP=NP+1  
      COEFN = A(NOEL,NP)
      NP=NP+1 
      COEFB = A(NOEL,NP)
      NP=NP+1 
      COEFG = A(NOEL,NP)

      ACN = ACENT * RAD      
      IF(NRES.GT.0) 
     >     WRITE(NRES,100)AT,ACENT,RM,HNORM,COEFN,COEFB,COEFG
C     >     WRITE(NRES,100)AT,ACENT,RM+DRM,HNORM,COEFN,COEFB,COEFG
  100 FORMAT(20X,'Dipole  magnet',//,1P,
     1 11X,'ANGLES : A.TOTAL =',G12.4,' degrees',5X,'A.CENTRAL =',
     2 G12.4,' degrees',/,11X,'RM =',G12.4,' cm',/,
C     2 G12.4,' degrees',/,11X,'RM+DRM =',G12.4,' cm',/,
     5 11X,'HNORM =',G12.4,' kGauss',5X,'COEF.N =',G12.4,5X,'COEF.B =',
     6 G12.4,5X,'COEF.G=',G12.4)

      NP=NP+1 
      LAMBDE = A(NOEL,NP)
      NP=NP+1 
      QSIE = A(NOEL,NP)
      NP=NP+1 
      NCOEFE = NINT(A(NOEL,NP))
      DO 227 I=1,6
        NP=NP+1 
 227    CE(KMAG,I) = A(NOEL,NP)
      NP=NP+1 
      SHIFTE = A(NOEL,NP)
      SHIFTE=0.D0
 
      NP=NP+1 
      UMEGA = A(NOEL,NP)
      NP=NP+1 
      THETA = A(NOEL,NP)
      NP=NP+1 
      R1    = A(NOEL,NP)
      NP=NP+1 
      U1    = A(NOEL,NP)
      NP=NP+1 
      U2    = A(NOEL,NP)
      NP=NP+1 
      R2    = A(NOEL,NP)
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,106) 'Entrance  EFB ',LAMBDE,QSIE
106     FORMAT (/,5X,A14,/,10X,
     >  'Fringe  field  : LAMBDA =',F7.2,' CM     QSI=',F7.2)
        WRITE(NRES,127) NCOEFE,(CE(KMAG,I),I=1,6),SHIFTE
127     FORMAT (10X,' COEFFICIENTS :',I3,6F10.5
     2  ,/,10X,' Shift  of  EFB  = ',G12.4,' CM',/)
        WRITE(NRES,103) UMEGA,THETA,R1,U1,U2,R2
C 103    FORMAT(10X,7HOMEGA =,F7.2,5X,17HANGLE  DE  FACE =,F7.2,/ ,
 103    FORMAT(10X,'OMEGA =',F7.2,' deg.',
     >   5X,'Wedge  angle  =',F7.2,' deg.',/,
     1   11X,'Radius 1 =',1P,G10.2,' CM',/ ,
     2   11X,'Straight  segment 1 =',G10.2,' CM',/ ,
     3   11X,'Straight  segment 2 =',G10.2,' CM',/ ,
     4   11X,'Radius 2 =',G10.2,' CM')
        IF(R1*R2 .EQ. 0.D0) WRITE(NRES,123)
 123    FORMAT(10('*'),' ATTENTION |    R1 OU R2 = 0 |',10('*'))
      ENDIF

C Exit Fringe Field
      NP=NP+1 
      LAMBDS = A(NOEL,NP)
      NP=NP+1 
      QSIS   = A(NOEL,NP)
      NP=NP+1 
      NCOEFS = NINT(A(NOEL,NP))
      DO 228 I=1,6
        NP=NP+1 
 228    CS(KMAG,I) = A(NOEL,NP)
      NP=NP+1 
      SHIFTS = A(NOEL,NP)
      SHIFTS=0D0

      NP=NP+1 
      UMEGAS = A(NOEL,NP)
      NP=NP+1 
      THETAS = A(NOEL,NP)
      NP=NP+1 
      R1S    = A(NOEL,NP)
      NP=NP+1 
      U1S    = A(NOEL,NP)
      NP=NP+1 
      U2S    = A(NOEL,NP)
      NP=NP+1 
      R2S    = A(NOEL,NP)
      IF(NRES.GT.0) THEN
        WRITE(NRES,106) 'Exit  EFB     ',LAMBDS,QSIS
        WRITE(NRES,127) NCOEFS,(CS(KMAG,I),I=1,6),SHIFTS
        WRITE(NRES,103) UMEGAS,THETAS,R1S,U1S,U2S,R2S
        IF(R1S*R2S .EQ. 0.D0) WRITE(NRES,123)
      ENDIF
 
      NP=NP+1 
      LAMBD3 = A(NOEL,NP)
      NP=NP+1 
      QSI3   = A(NOEL,NP)
      NP=NP+1 
      NCOEF3 = NINT(A(NOEL,NP))
      DO 229 I=1,6
        NP=NP+1 
 229    C3(KMAG,I) = A(NOEL,NP)
      NP=NP+1 
      SHIFT3 = A(NOEL,NP)
      SHIFT3=0D0
 
       NP=NP+1 
       UMEGA3 = A(NOEL,NP)
       NP=NP+1 
       THETA3 = A(NOEL,NP)
       NP=NP+1 
       R13    = A(NOEL,NP)
       NP=NP+1 
       U13    = A(NOEL,NP)
       NP=NP+1 
       U23    = A(NOEL,NP)
       NP=NP+1 
       R23    = A(NOEL,NP)
       NP=NP+1 
       RM3    = A(NOEL,NP)
 
       IF(NRES.GT.0) THEN
         WRITE(NRES,106) 'Lateral  EFB  ',LAMBD3,QSI3
         WRITE(NRES,127) NCOEF3,(C3(KMAG,I),I=1,6),SHIFT3
         WRITE(NRES,113) RM3
113      FORMAT(20X,' Face centred on direction ACENT+OMEGA, A ',
     >   G12.4,' CM'/)
         WRITE(NRES,103) UMEGA3,THETA3,R13,U13,U23,R23
         IF(R13*R23 .EQ. 0.D0) WRITE(NRES,123)
       ENDIF
       
      NP=NP+1 
      IRD = NINT(A(NOEL,NP))
      IF    (IRD.EQ.2) THEN 
        NN=3
      ELSEIF(IRD.EQ.25) THEN
        NN=5
      ELSEIF(IRD.EQ.4) THEN
        NN=5
      ELSE
        STOP ' *** ERROR - SBR DIPI, WRONG VALUE IRD'
      ENDIF
      NP=NP+1 
      RESOL=A(NOEL,NP)

      AE=0.D0
      AS=0.D0
      AT = AT * RAD
      XI = 0.D0
      XF = AT

      DSREF = RM * ((UMEGA-UMEGAS)*RAD +  TAN(ACN - UMEGA* RAD) + 
     >         TAN(AT - ACN + UMEGAS* RAD) )    

      IF(NRES.GT.0) THEN
        WRITE(NRES,119) IRD
 119    FORMAT(/,20X,' Interpolation  option :',I2)
        IF(IRD .EQ. 2) THEN
          WRITE(NRES,121) ' 9', RESOL
 121      FORMAT(20X,A2,
     >     '-point  interpolation, size of flying mesh :   STEP /',F7.1)
        ELSE
C--------- IRD=4 OU 25
          WRITE(NRES,121) '25', RESOL
        ENDIF
      ENDIF

      RETURN

      ENTRY DIPF(TTA,RO,
     >                  AAA,RRR,DA,DR,FTAB)
      CALL INTEG5(
     >            STEP)
      DRO = STEP/RESOL
      DTTA = DRO/RM
      DR = DRO
      DA = DTTA

C Entrance EFB
C     ** POUR LES BESOINS DE LA SIMPLE PRECISION :
      IF(R1*R1 .GE. 1.D6) U1 = -1.D6
      IF(R2*R2 .GE. 1.D6) U2 =  1.D6
 
      UMEG = UMEGA * RAD
      TETA = THETA * RAD
      UT = UMEG - TETA
      SINO = SIN(UT)
      COSO = COS(UT)
      TANO = SINO / COSO
 
C  CALCUL DES PARAMETRES DE FACE  :
C  AXE X DU REFERENTIEL  = PARALLELE  DIRCTN ACENT, SENS DES RAYONS CROISSANTS
C  AXE Y DU REFERENTIEL  = ORTHOGONAL DIRCTN ACENT, DIRIGE VERS FACE D'ENTREE
C  ORIGINE DU REFERENTIEL= A L'INTERSECTION DE RM ET ACENT
 
C  PROJECTIONS DU RAYON MOYEN SUR LES AXES X,Y
      XB = RM * ( COS(UMEG) - 1.D0)
      YB = RM * SIN(UMEG)
C  COORDONNEES DE L'EXTREMITE A DE LA PARTIE LINEAIRE DE LONG. U2
      XA = U2 * COSO + XB
      YA = YB + U2 * SINO
C  COORDONNEES DE L'EXTREMITE C DE LA PARTIE LINEAIRE DE LONG. U1
      XC = U1 * COSO + XB
      YC = YB + U1 * SINO
C  COORDONNEES DU CENTRE DE COURBURE DE RAYON R1
      XD = XC + R1 * SINO
      YD = YC - R1 * COSO
C  COORDONNEES DU CENTRE DE COURBURE DE RAYON R2
      XE = XA + R2 * SINO
      YE = YA - R2 * COSO
      SIN2 = SINO**2
      COS2 = COSO**2
      SICO = SINO * COSO
 
C Exit EFB
C     ** POUR LES BESOINS DE LA SIMPLE PRECISION :
      IF(R1S*R1S .GE. 1.D10) U1S = -1.D6
      IF(R2S*R2S .GE. 1.D10) U2S =  1.D6
 
      UMEGS = UMEGAS * RAD
      TETAS = THETAS * RAD
      UTS= UMEGS- TETAS
      SINOS= SIN(UTS)
      COSOS= COS(UTS)
      TANOS= SINOS/ COSOS
      XBS = RM * ( COS(UMEGS) - 1.D0)
      YBS = RM * SIN(UMEGS)
      XAS = U2S * COSOS + XBS
      YAS = YBS + U2S * SINOS
      XCS = U1S * COSOS + XBS
      YCS = YBS + U1S * SINOS
      XDS = XCS + R1S * SINOS
      YDS = YCS - R1S * COSOS
      XES = XAS + R2S * SINOS
      YES = YAS - R2S * COSOS
      SIN2S = SINOS**2
      COS2S = COSOS**2
      SICOS = SINOS * COSOS
 
C Third EFB
      NBFACE=2
      IF(LAMBD3 .NE. 0.D0) NBFACE=3

      IF(NBFACE .EQ. 3) THEN
        UMEG3 = UMEGA3 * RAD
        TETA3 = THETA3 * RAD
        UT3= UMEG3- TETA3
        IF(ABS(UT3) .GE. .5D0*PI) THEN
          IF(NRES.GT.0) WRITE(NRES,139)
139       FORMAT(/,20X,' FACE 3, ATTENTION : ABS(UT3) > 90 DEG. =>',
     >    ' VALEUR AMBIGUE POUR LE CALCUL DE LA CARTE DE Champ')
        ENDIF
 
        SINO3= SIN(UT3)
        COSO3= COS(UT3)
        TANO3= SINO3/ COSO3
        XB3 = RM3* COS(UMEG3) - RM
        YB3 = RM3* SIN(UMEG3)
        XA3 = U23 * COSO3 + XB3
        YA3 = YB3 + U23 * SINO3
        XC3 = U13 * COSO3 + XB3
        YC3 = YB3 + U13 * SINO3
        XD3 = XC3 + R13 * SINO3
        YD3 = YC3 - R13 * COSO3
        XE3 = XA3 + R23 * SINO3
        YE3 = YA3 - R23 * COSO3
        SIN23 = SINO3**2
        COS23 = COSO3**2
        SICO3 = SINO3 * COSO3
      ENDIF
C---------- endif third EFB

C----- CALCUL LE Champ en X, Y 
C      COORDONNEES DU POINT COURANT
      RRR = 0.D0
      DO  1  JRO = 1,NN
        ROJ = RO + DRO * DBLE(NN-JRO-INT(NN/2))

        LTXI = TTA-DTTA .LT. XI
        GTXF = TTA+DTTA .GT. XF
        OKTTA = .NOT.(LTXI .OR. GTXF)
C        write(*,*) ' dipi  oktta, ltxi, gtxf ', oktta, ltxi, gtxf
        IF(OKTTA) THEN 
          AAA = 0.D0
        ELSEIF(LTXI) THEN 
          AAA=-DTTA 
          DA = 2.D0 * DTTA
          IF(IRD .NE. 2) THEN
            AAA = AAA -DTTA
            DA = 4.D0 * DTTA
          ENDIF
        ELSEIF(GTXF) THEN
          AAA= DTTA 
          DA = 2.D0 * DTTA
          IF(IRD .NE. 2) THEN
            AAA = AAA +DTTA
            DA = 4.D0 * DTTA
          ENDIF          
        ENDIF

        DO  1  ITTA = 1,NN
          TTAI = TTA - DTTA * DBLE(NN-ITTA-INT(NN/2))
          IF(OKTTA) THEN 
          ELSEIF(LTXI) THEN 
            TTAI = TTAI + DTTA
            IF(IRD .NE. 2) TTAI = TTAI + DTTA
          ELSEIF(GTXF) THEN
            TTAI = TTAI - DTTA
            IF(IRD .NE. 2) TTAI = TTAI - DTTA
          ENDIF
          
          ZETA = ACN - TTAI 

          X = ROJ * COS(ZETA) - RM
          Y = ROJ * SIN(ZETA)

C         ... AX (CX) = COTE X DE LA PROJECTION DE A (C) SUR LA DROITE (MX)
C         (M=POINT COURANT), ORTHOGONALEMENT A LA DROITE (ABC)
          AX = ((YA - Y) * TANO ) + XA
          CX = ((YC - Y) * TANO ) + XC
          AXS = ((YAS - Y) * TANOS ) + XAS
          CXS = ((YCS - Y) * TANOS ) + XCS
          IF(NBFACE .EQ. 3) THEN
            AX3 = ((YA3 - Y) * TANO3 ) + XA3
            CX3 = ((YC3 - Y) * TANO3 ) + XC3
          ENDIF
C
C POSITION DU POINT (X,Y) / FACE ENTREE
          IF     ( X.LT.CX .AND. X.LT.AX )  THEN
C           REGION DE COURBURE R1
            RR = R1
            XO = XD
            YO = YD
            R = ABS(RR)
            D = R - SQRT((X - XO)**2 + (Y - YO)**2 )
            D =( -D * RR/R + SHIFTE )
          ELSE IF( X.GE.CX .AND. X.LE.AX )  THEN
C            REGION LINEAIRE
            XO = SICO * (Y - YB) + XB * SIN2 + X * COS2
            YO = SICO * (X - XB) + YB * COS2 + Y * SIN2
            YL = YB + ((X - XB) * TANO )
C           ... YL = COTE Y DE LA PROJECTION DE B SUR LA DROITE (MY)
C           (M=POINT COURANT), PARALLELEMENT A LA DROITE (ABC)
            D = SQRT((X - XO)**2 + (Y - YO)**2 )
            IF( Y .LE. YL .OR. D .LE. 1.D-6 ) D = -D
            D=( D + SHIFTE )
          ELSE IF( X.GT.CX .AND. X.GT.AX )  THEN
C           REGION DE COURBURE R2
            RR = R2
            XO = XE
            YO = YE
            R = ABS(RR)
            D = R - SQRT((X - XO)**2 + (Y - YO)**2 )
            D=( -D * RR/R + SHIFTE )
          ELSE
C           ERREUR  DE  DONNEES  FACE  ENTREE
            IF(NRES.GT.0) WRITE(NRES,104)
 104        FORMAT(/,5X,10('*'),' ERREUR PARAMETRES FACE ENTREE',/)
            IF(NRES .GT. 0) WRITE(NRES,134) X,AX,CX,XA,YA,XC,YC
            GOTO  99
          ENDIF
 
          IF(LAMBDE .EQ. 0.D0) THEN
            IF(D.LE.0.D0) THEN
              FE=1.D0
            ELSE
              FE=0.D0
            ENDIF
          ELSE
            D = D/LAMBDE
            P=CE(KMAG,1)+(CE(KMAG,2)+(CE(KMAG,3)+(CE(KMAG,4)+
     >                       (CE(KMAG,5)+CE(KMAG,6)*D)*D)*D)*D)*D
            IF    (P .GE.  12.D0) THEN
              FE = 0.D0
            ELSEIF(P .LE. -12.D0) THEN
              FE = 1.D0
            ELSE
              FE = 1.D0/(1.D0+EXP(P))
            ENDIF
          ENDIF
 
C POSITION DU POINT (X,Y) / FACE SORTIE
          IF     ( X.LT.CXS .AND. X.LT.AXS )  THEN
C           ... REGION DE COURBURE R1
            RRS = R1S
            XO = XDS
            YO = YDS
            R = ABS(RRS)
            D = R - SQRT((X - XO)**2 + (Y - YO)**2 )
            D=( D * RRS/R + SHIFTS )
          ELSE IF( X.GE.CXS .AND. X.LE.AXS )  THEN
C           ... REGION LINEAIRE
            XO = SICOS * (Y - YBS) + XBS * SIN2S + X * COS2S
            YO = SICOS * (X - XBS) + YBS * COS2S + Y * SIN2S
            YL = YBS + ((X - XBS) * TANOS )
            D = SQRT((X - XO)**2 + (Y - YO)**2 )
            IF( Y .GE. YL .OR. D .LE.1.D-6 ) D = -D
            D=( D + SHIFTS )
          ELSE IF( X.GT.CXS .AND. X.GT.AXS )  THEN
C           ... REGION DE COURBURE R2
            RRS = R2S
            XO = XES
            YO = YES
            R = ABS(RRS)
            D = R - SQRT((X - XO)**2 + (Y - YO)**2 )
            D=( D * RRS/R + SHIFTS )
          ELSE
            IF(NRES.GT.0) WRITE(NRES,114)
114         FORMAT(/,5X,10('*'),' ERREUR PARAMETRES FACE SORTIE',/)
            IF(NRES .GT. 0) WRITE(NRES,134) X,AXS,CXS,XAS,YAS,XCS,YCS
            GOTO  99
          ENDIF
 
          IF(LAMBDS .EQ. 0.D0) THEN
            IF(D.LE.0.D0) THEN
              FS=1.D0
            ELSE
              FS=0.D0
            ENDIF
          ELSE
            D = D/LAMBDS
            P=CS(KMAG,1)+(CS(KMAG,2)+(CS(KMAG,3)+(CS(KMAG,4)+
     >                        (CS(KMAG,5)+CS(KMAG,6)*D)*D)*D)*D)*D
            IF    (P .GE.  14.D0) THEN
              FS = 0.D0
            ELSEIF(P .LE. -14.D0) THEN
              FS = 1.D0
            ELSE
              FS = 1.D0/(1.D0+EXP(P))
            ENDIF
          ENDIF
 
C POSITION DU POINT (X,Y) / FACE EXTERNE (SEULEMENT SI OPTION 3 FACES)
          IF(NBFACE .EQ. 3) THEN
            IF     ( X.LT.CX3 .AND. X.LT.AX3 )  THEN
C             ... REGION DE COURBURE R1
              RR3 = R13
              XO = XD3
              YO = YD3
              R = ABS(RR3)
              D = R - SQRT((X - XO)**2 + (Y - YO)**2 )
              D=( D * RR3/R + SHIFT3 )
            ELSE IF( X.GE.CX3 .AND. X.LT.AX3 )  THEN
C             ... REGION LINEAIRE
              XO = SICO3 * (Y - YB3) + XB3 * SIN23 + X * COS23
              YO = SICO3 * (X - XB3) + YB3 * COS23 + Y * SIN23
              YL = YB3 + ((X - XB3) * TANO3 )
              D = SQRT((X - XO)**2 + (Y - YO)**2 )
              IF( Y .GE. YL .OR. D .LE. 1.D-6 ) D = -D
              D=( D + SHIFT3 )
            ELSE IF( X.GE.CX3 .AND. X.GE.AX3 )  THEN
C             ... REGION DE COURBURE R2
              RR3 = R23
              XO = XE3
              YO = YE3
              R = ABS(RR3)
              D = R - SQRT((X - XO)**2 + (Y - YO)**2 )
              D=( D * RR3/R + SHIFT3 )
            ELSE
              IF(NRES.GT.0) THEN
                WRITE(NRES,124)
124             FORMAT(/,5X,9('*'),'ERREUR PARAMETRES FACE LATERALE',/)
                WRITE(NRES,134) X,AX3,CX3,XA3,YA3,XC3,YC3
134             FORMAT(/,5X,7G12.4)
              ENDIF
              GOTO  99
            ENDIF
 
            IF(LAMBD3 .EQ. 0.D0) THEN
              IF(D.LE.0.D0) THEN
                F3=1.D0
              ELSE
                F3=0.D0
              ENDIF
            ELSE
              D = D/LAMBD3
              P=C3(KMAG,1)+(C3(KMAG,2)+(C3(KMAG,3)+(C3(KMAG,4)+
     >                           (C3(KMAG,5)+C3(KMAG,6)*D)*D)*D)*D)*D
              IF    (P .GE.  14.D0) THEN
                F3 = 0.D0
              ELSEIF(P .LE. -14.D0) THEN
                F3 = 1.D0
              ELSE
                F3 = 1.D0/(1.D0+EXP(P))
              ENDIF
            ENDIF
          ENDIF
 
          F = FE * FS
          IF(NBFACE .EQ. 3) F = F * F3
 
         ROI= (ROJ-RM)/RM
         FTAB(ITTA,JRO)=F*HNORM*(1.D0+(COEFN+(COEFB+COEFG*ROI)*ROI)*ROI)

    1 CONTINUE
C--------- end loop on   NN x tta  &  NN x ro

      RETURN
 99   STOP
      END 
