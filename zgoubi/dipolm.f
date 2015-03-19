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
      SUBROUTINE DIPOLM(SCAL,
     >                          BMIN,BMAX,BNORM,
     >                          XBMI,YBMI,ZBMI,XBMA,YBMA,ZBMA)
      USE dynhc
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     Build a simulated 2D magnetic field map in polar coordinates mesh.
C     Improved version of SBR AIMANT.
C-----------------------------------------------------------------------
      DOUBLE PRECISION LAMBDS,LAMBDE,LAMBD3
      INCLUDE 'PARIZ.H'
      INCLUDE "XYZHC.H"
      INCLUDE "C.AIM_2.H"     ! COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CARSH.H"     ! COMMON/CARSH/ ATS,RMS
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH 
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.DROITE_2.H"     ! COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
 
      DIMENSION CI(6),CS(6),C3(6)
C      DATA IMAP / 1 /

      CALL KSMAP5(
     >            IMAP)
      if(imap.le.0 .or. imap .gt. mmap)
     >call endjob('Pgm dipolm. Wrong value IMAP = ',imap)

C----- NBFACE=2(3) :DIPOLE LIMITE PAR (2)3 FACES.
      NBFACE = NINT(A(NOEL,1))
 
C  IXMA/JYMA = NBRE DE PAS ANGULAIRE/RADIAL DE LA CARTE DE Champ.
C  HNORM=Champ MAX DANS LE DIPOLE.
C  COEFN=N=INDICE DE Champ, B=N', G=N''.
C  AT=ANGLE TOTAL DE LA CARTE DE Champ, ACENT='ANGLE AU CENTRE',
C    RM,MIN,MAX=RAYONS MOYEN,MIN,MAX DE LACARTE DE Champ.
 
      IXMA = NINT(A(NOEL,4))
      IF(IXMA.GT.MXX) 
     >   CALL ENDJOB('X-dim of map is too large, max is ',MXX)
      JYMA = NINT(A(NOEL,5))
      IF(JYMA.GT.MXY ) 
     >   CALL ENDJOB('Y-dim of map is too large, max is ',MXY)

      HNORM = A(NOEL,6)*SCAL
      COEFN = A(NOEL,7)
      COEFB = A(NOEL,8)
      COEFG = A(NOEL,9)

      AT    = A(NOEL,10)
      ACENT = A(NOEL,11)
      RM    = A(NOEL,12)
      RMIN  = A(NOEL,13)
      RMAX  = A(NOEL,14)
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,100) AT,ACENT,RM,RMIN,RMAX,HNORM,COEFN,COEFB,COEFG
        IF(NBFACE .EQ. 3 .AND. (RMIN.GT.RM.OR.RMAX.LT.RM))
     >  WRITE(NRES,101)
101     FORMAT(/,10X,' ATTENTION: VALEURS DES RAYONS INCOMPATIBLES||'/)
        WRITE(NRES,105) IXMA,JYMA
105     FORMAT(20X,'NBRE PAS ANGULAIRE=',I5,'  NBRE PAS RAYON=',I5)
      ENDIF
 
      AE = 0.D0
      AS = 0.D0
      ATS = AT / (IXMA - 1)
      RMS = (RMAX - RMIN) / (JYMA - 1)
      DO  25  I = 1,JYMA
         YH(I) = RMIN + (I - 1) * RMS
   25 CONTINUE
      NN = 0
      ACN = ACENT * RAD
      ATSRM=ATS*RAD*RM
      IF(NRES.GT.0) WRITE(NRES,120)ATS,ATSRM,RMS
120   FORMAT(19X,' PAS ANGULAIRE =',F8.4,' DEG. OU ',F8.4
     1,' CM AU RAYON MOYEN',/,19X,' PAS EN RAYON  =',F8.4,' CM')
 
C LE MODELE DE Champ DE FUITE DE CHACUNE DES 3 FACES EST CELUI DE ENGE
C Champ DE FUITE ENTREE
      LAMBDE = A(NOEL,15)
      DO 227 I=1,6
 227    CI(I) = A(NOEL,17+I)
      SHIFTE = A(NOEL,24)
      SHIFTE=0.D0
 
      UMEGA = A(NOEL,25)
      THETA = A(NOEL,26)
      R1    = A(NOEL,27)
      U1    = A(NOEL,28)
      U2    = A(NOEL,29)
      R2    = A(NOEL,30)
 
C     ** POUR LES BESOINS DE LA SIMPLE PRECISION :
      IF(R1*R1 .GE. 1.D6) U1 = -1.D6
      IF(R2*R2 .GE. 1.D6) U2 =  1.D6
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,106) LAMBDE
106     FORMAT (/,5X,'Entrance  face',/,10X,
C106     FORMAT (/,5X,'FACE  D''ENTREE',/,10X,
     >  'Champ FUITE : LAMBDA =',F7.2,' CM')
        WRITE(NRES,127) (CI(I),I=1,6),SHIFTE
127     FORMAT (10X,' COEFFICIENTS :',6F10.5
     2  ,/,10X,' DECALAGE  FACE  MAGNETIQUE = ',G12.4,' CM',/)
        WRITE(NRES,103) UMEGA,THETA,R1,U1,U2,R2
        IF(R1*R2 .EQ. 0.D0) WRITE(NRES,123)
 123    FORMAT(10('*'),' ATTENTION |    R1 OU R2 = 0 |',10('*'))
      ENDIF
 
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
 
C Champ DE FUITE SORTIE
      LAMBDS = A(NOEL,31)
      DO 228 I=1,6
 228    CS(I) = A(NOEL,33+I)
      SHIFTS = A(NOEL,40)
      SHIFTS=0D0
 
      UMEGAS = A(NOEL,41)
      THETAS = A(NOEL,42)
      R1S    = A(NOEL,43)
      U1S    = A(NOEL,44)
      U2S    = A(NOEL,45)
      R2S    = A(NOEL,46)
 
C     ** POUR LES BESOINS DE LA SIMPLE PRECISION :
      IF(R1S*R1S .GE. 1.D10) U1S = -1.D6
      IF(R2S*R2S .GE. 1.D10) U2S =  1.D6
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,116) LAMBDS
116     FORMAT(/,5X,'Exit  face',/,10X,
C116     FORMAT(/,5X,'FACE  DE  SORTIE',/,10X,
     >  'Champ FUITE : LAMBDA =',F7.2,' CM')
        WRITE(NRES,127) (CS(I),I=1,6),SHIFTS
        WRITE(NRES,103) UMEGAS,THETAS,R1S,U1S,U2S,R2S
        IF(R1S*R2S .EQ. 0.D0) WRITE(NRES,123)
      ENDIF
 
      UMEGS = UMEGAS * RAD
      TETAS = THETAS * RAD
      UTS= UMEGS- TETAS
      SINOS= SIN(UTS)
      COSOS= COS(UTS)
      TANOS= SINOS/ COSOS
      XBS = RM * ( COS(UMEGS) - 1D0 )
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
 
C Champ DE FUITE FACE EXTERNE (SEULEMENT SI OPTION 3 FACES)
      IF(NBFACE .EQ. 3) THEN
       LAMBD3 = A(NOEL,47)
       DO 229 I=1,6
 229     C3(I) = A(NOEL,49+I)
       SHIFT3 = A(NOEL,56)
       SHIFT3=0D0
 
       UMEGA3 = A(NOEL,57)
       THETA3 = A(NOEL,58)
       R13    = A(NOEL,59)
       U13    = A(NOEL,60)
       U23    = A(NOEL,61)
       R23    = A(NOEL,62)
       RM3    = A(NOEL,63)
 
       IF(NRES.GT.0) THEN
         WRITE(NRES,136) LAMBD3
136      FORMAT (/,5X,'FACE  LATERALE ',/,10X,
     >   'Champ FUITE : LAMBDA =',F7.2,' CM')
         WRITE(NRES,127) (C3(I),I=1,6),SHIFT3
         WRITE(NRES,113) RM3
113      FORMAT(20X,' FACE CENTREE SUR LA DIRECTION ACENT+OMEGA, A ',
     >   G12.4,' CM'/)
         WRITE(NRES,103) UMEGA3,THETA3,R13,U13,U23,R23
         IF(R13*R23 .EQ. 0.D0) WRITE(NRES,123)
       ENDIF
 
       UMEG3 = UMEGA3 * RAD
       TETA3 = THETA3 * RAD
       UT3= UMEG3- TETA3
       IF(ABS(UT3) .GE. .5D0*PI) THEN
         IF(NRES.GT.0) WRITE(NRES,139)
139      FORMAT(/,20X,' FACE 3, ATTENTION : ABS(UT3) > 90 DEG. =>',
     >   ' VALEUR AMBIGUE POUR LE CALCUL DE LA CARTE DE Champ')
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
 
C  CALCUL LA CARTE DE Champ, POUR IXMA ANGLES ET JYMA RAYONS
 
C     ****DEBUT DE BOUCLES SUR ANGLES ET RAYONS
      BMIN=1.D10
      BMAX=-1.D10
      ZBMI = 0D0
      ZBMA = 0D0

      DO  1  J = 1,IXMA
        NN = NN + 1
        XH(NN) = (NN - 1) * ATS * RAD
        ZETA = ACN - XH(NN)
        DO  2  I = 1,JYMA
 
C         ... COORDONNEES DU POINT COURANT
          X = YH(I) * COS(ZETA) - RM
          Y = YH(I) * SIN(ZETA)

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
            GOTO  15
          ENDIF
 
          IF(LAMBDE .EQ. 0.D0) THEN
            IF(D.LE.0.D0) THEN
              FE=1D0
            ELSE
              FE=0D0
            ENDIF
          ELSE
            D = D/LAMBDE
            P=CI(1)+(CI(2)+(CI(3)+(CI(4)+(CI(5)+CI(6)*D)*D)*D)*D)*D
            IF    (P .GE.  12.D0) THEN
              FE = 0D0
            ELSEIF(P .LE. -12.D0) THEN
              FE = 1D0
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
            GOTO  15
          ENDIF
 
          IF(LAMBDS .EQ. 0.D0) THEN
            IF(D.LE.0.D0) THEN
              FS=1D0
            ELSE
              FS=0D0
            ENDIF
          ELSE
            D = D/LAMBDS
            P=CS(1)+(CS(2)+(CS(3)+(CS(4)+(CS(5)+CS(6)*D)*D)*D)*D)*D
            IF    (P .GE.  14.D0) THEN
              FS = 0D0
            ELSEIF(P .LE. -14.D0) THEN
              FS = 1D0
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
              GOTO  15
            ENDIF
 
            IF(LAMBD3 .EQ. 0.D0) THEN
              IF(D.LE.0.D0) THEN
                F3=1D0
              ELSE
                F3=0D0
              ENDIF
            ELSE
              D = D/LAMBD3
              P=C3(1)+(C3(2)+(C3(3)+(C3(4)+(C3(5)+C3(6)*D)*D)*D)*D)*D
              IF    (P .GE.  14.D0) THEN
                F3 = 0D0
              ELSEIF(P .LE. -14.D0) THEN
                F3 = 1D0
              ELSE
                F3 = 1.D0/(1.D0+EXP(P))
              ENDIF
            ENDIF
          ENDIF
 
          F = FE * FS
          IF(NBFACE .EQ. 3) F = F * F3
 
C         CE 'POINT DE MESURE' PERMET DE SUIVRE LE
C         PARCOURS A TRAVERS LES FACES, LE LONG D'UN RAYON.
C         IF(I .EQ. JYMA/2)
C     X   WRITE(NRES,8001) FE,FS,F3,F,X,Y,J,I
C 8001    FORMAT(' .....FE,FS,F3,F  ,X,Y, ANG,RAY :',6G12.4,2I5)
C
          YHI= (YH(I)-RM)/RM
          BFLD=F*HNORM*(1.D0+(COEFN+(COEFB+COEFG*YHI)*YHI)*YHI)
          IF    (BFLD .GT. BMAX) THEN
            BMAX = BFLD
            XBMA = XH(NN)
            YBMA = YH(I)
          ELSEIF(BFLD .LT. BMIN) THEN
            BMIN = BFLD
            XBMI = XH(NN)
            YBMI = YH(I)
          ENDIF
          HC(ID,NN,I,1,IMAP) = BFLD
C
C       ... FIN DE BOUCLES SUR ANGLES ET RAYONS
    2   CONTINUE
    1 CONTINUE

      BNORM=HNORM

C     ****SIMULE DES PERTURBATIONS DE LA CARTE DE Champ
      IF(NBFACE .EQ. 2) THEN
        ND = 57
      ELSEIF(NBFACE .EQ. 3) THEN
        ND = 64
      ENDIF
      NBSHIM = NINT(A(NOEL,ND))
      IF    (NBSHIM .NE. 0) THEN
         IF(NBSHIM.GT.0) THEN
C           ****INTRODUCTION D'ILOTS OU SHIMS
            IF(NRES.GT.0) 
     >      WRITE(NRES,*) ' NO SHIMS AVAILABLE: TO BE IMPLEMENTED'
C            CALL CARSHI(NBSHIM)
         ELSEIF(NBSHIM .EQ. -1) THEN
C           ****PERTURBATION LINEAIRE EN ANGLE
            ZETA = A(NOEL,ND+1)
            DBSB = A(NOEL,ND+2)
            IF(NRES.GT.0) WRITE(NRES,109) ZETA,DBSB
109         FORMAT(/,20X,' PERTURBATION DE LA CARTE DE Champ, LINEAIRE'
     >      ,'EN ANGLE ET',/,20X,' CENTREE SUR THETA =',F7.2
     >      ,' DEG.,  AVEC DB/B0(TOTALE) = ',F10.6)
            J0=INT(ZETA/ATS)
            DZETA=ZETA-J0*ATS
            DBSB=DBSB /IXMA
            P=1.D0-DBSB*DZETA/ATS
            DO 3 J=1,IXMA
              DO 3 I=1,JYMA
                HC(ID,J ,I,1 ,IMAP)=HC(ID,J ,I,1 ,IMAP)*(P+DBSB*(J-J0))
3           CONTINUE
         ELSEIF(NBSHIM .EQ. -2) THEN
C           ****PERTURBATION LINEAIRE EN RAYON
            RR   = A(NOEL,ND+1)
            DBSB = A(NOEL,ND+2)
            IF(NRES.GT.0) WRITE(NRES,110) RR  ,DBSB
110         FORMAT(/,20X,' PERTURBATION DE LA CARTE DE Champ, LINEAIRE'
     >      ,'EN RAYON ET',/,20X,' CENTREE SUR R =',F7.2
     >      ,' CM,   AVEC DB/B0(TOTALE) = ',F10.6)
            J0=INT(RR  /RMS)
            DRR  =RR  -J0*RMS
            DBSB=DBSB /JYMA
            P=1.D0-DBSB*DRR  /RMS
            DO 4 J=1,IXMA
              DO 4 I=1,JYMA
                HC(ID,J ,I,1 ,IMAP)=HC(ID,J ,I,1 ,IMAP)*(P+DBSB*(I-J0))
4           CONTINUE
         ENDIF
      ENDIF
 
      BAMP = 1.D10
      CALL MAPLI1(BAMP)
 
      IRD = NINT(A(NOEL,93))
 
      AT = AT * RAD
      XI = XH(1)
      XF = XH(NN)
 
15    CONTINUE
 
  100 FORMAT(20X,'AIMANT  PRINCIPAL',//,
     1 11X,'ANGLES : A.TOTAL =',F6.2,' degrees',5X,'A.CENTRAL =',
     2 F6.2,' degrees',/ ,11X,'RM =',F7.2,' cm',5X,'RMIN =',F7.2,' cm',
     3 5X,'RMAX =',F7.2,' cm',/,
     5 11X,'HNORM =',F8.4,' kGauss',5X,'COEF.N =',F9.5,5X,'COEF.B =',
     6F9.5,5X,'COEF.G=',F9.5)
C  103 FORMAT(10X,7HOMEGA =,F7.2,5X,17HANGLE  DE  FACE =,F7.2,/ ,
 103  FORMAT(10X,'OMEGA =',F7.2,' deg.',
     > 5X,'Wedge  angle  =',F7.2,' deg.',/,
     1 11X,'Radius 1 =',1P,G10.2,' cm',/ ,
     2 11X,'Straight  segment 1 =',G10.2,' cm',/ ,
     3 11X,'Straight  segment 2 =',G10.2,' cm',/ ,
     4 11X,'Radius 2 =',G10.2,' cm')
      RETURN
      END
