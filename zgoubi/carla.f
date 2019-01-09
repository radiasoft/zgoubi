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
      SUBROUTINE CARLA(SCAL,
     >                          BMIN,BMAX,BNORM,
     >                          ABMI,RBMI,ZBMI,ABMA,RBMA,ZBMA)
      USE dynhc
      use xyzhc_interface, only : XH, YH, IXMA, JYMA
      use pariz_namelist_interface, only : ID, MXX, MXY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     *******
C     CONSTRUIT LA CARTE DE Champ D'UN DIPOLE EN MAILLAGE
C       POLAIRE. LE DIPOLE PRESENTE 2 (3) FACES: ENTREE, SORTIE
C       (FACE EXTERIEURE). CHAQUE FACE PEUT ETRE CONSTITUEE D'UNE
C       PARTIE RECTILIGNE ET DE 2 PARTIES COURBES.
C     CES 2 (3) FACES SONT DOTEES D'UN Champ DE FUITE DE TYPE SECOND
C       ORDRE, OU BIEN B=B0/(1+EXP(P(X))) OU P(X) EST UN POLYNOME
C       DE DEGRE AU PLUS EGAL A 6.
C     *******
      DOUBLE PRECISION LAMBDS,LAMBDE,LDRTS,LDRTE,LO
      DOUBLE PRECISION LAMBD3,LDRT3
      INCLUDE 'C.AIM_2.H'     ! COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CARSH.H"     ! COMMON/CARSH/ ATS,RMS
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"

      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.DROITE_2.H"     ! COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      INCLUDE "C.ORDRES.H"     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI

      DIMENSION CI(6),CS(6),C3(6)

      CALL KSMAP5(
     >            IMAP)

C  NBFACE=(2)3 :DIPOLE LIMITE PAR (2)3 FACES.
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


C
      IF(NRES.GT.0) THEN
        WRITE(NRES,100) AT,ACENT,RM,RMIN,RMAX,HNORM,COEFN,COEFB,COEFG
        IF(NBFACE .EQ. 3 .AND. (RMIN.GT.RM.OR.RMAX.LT.RM))
     >  WRITE(NRES,101)
101     FORMAT(/,10X,' ATTENTION : Incompatible  radius  values !!'/)
C101     FORMAT(/,10X,' ATTENTION: VALEURS DES RAYONS INCOMPATIBLES||'/)
        WRITE(NRES,105) IXMA,JYMA
105     FORMAT(20X,'Number  of  steps,  in  angle :',
     >                                    I5,'  in  radius :',I5)
C105     FORMAT(20X,'NBRE PAS ANGULAIRE=',I5,'  NBRE PAS RAYON=',I5)
      ENDIF

      AE = 0D0
      AS = 0D0
      ATS = AT / (IXMA - 1)
      RMS = (RMAX - RMIN) / (JYMA - 1)
      DO  25  I = 1,JYMA
         YH(I) = RMIN + (I - 1) * RMS
   25 CONTINUE
      NN = 0
      ACN = ACENT * RAD
      ATSRM=ATS*RAD*RM
      IF(NRES.GT.0) WRITE(NRES,120)ATS,ATSRM,RMS
120   FORMAT(19X,' Angular  step =',F8.4,' deg.,  i.e.,  ',F8.4
     1,' cm  at  mean  radius',/,19X,' Radial  step  =',F8.4,' cm')
C120   FORMAT(19X,' PAS ANGULAIRE =',F8.4,' DEG. OU ',F8.4
C     1,' CM AU RAYON MOYEN',/,19X,' PAS EN RAYON  =',F8.4,' CM')
C
C LE MODELE DE Champ DE FUITE DE CHACUNE DES 3 FACES
C   EST DETERMINE PAR QSIE/QSIS/QSI3:
C   -QSIE/QSIS/3.LT.0. : Champ DE FUITE ENTREE/SORTIE EN 1/(1+EXP(P))
C   ET INDICATEUR IQSIE/IQSIS/IQSI3=0,
C   -SINON Champ DE FUITE ORDRE 2
C   ET INDICATEUR IQSIE/IQSIS/IQSI3=1.
C LAMBDA EST LA DEMI-ETENDUE DU Champ DE FUITE, ET LES COEFFS SONT LES
C   COEFFICIENTS DU MODELE POLYNOMIAL (CES COEFFICIENTS SONT INOPERANTS
C   SI LE Champ DE FUITE ORDRE 2 EST CHOISI)
C
C Champ DE FUITE ENTREE
      LAMBDE = A(NOEL,15)
      QSIE   = A(NOEL,16)
      NCOEFE = NINT(A(NOEL,17))
      DO 227 I=1,6
 227    CI(I) = A(NOEL,17+I)
      SHIFTE = A(NOEL,24)
      SHIFTE=0D0

      NCFE1=NCOEFE+1
      IF(QSIE.LT.0.D0) THEN
         IF(NRES.GT.0) WRITE(NRES,107)
107      FORMAT(/2X,' Exponential  entrance  fringe  field'/)
C107      FORMAT(/2X,' Champ DE FUITE ENTREE EXPONENTIEL'/)
         QSIE=0D0
         IQSIE=0
      ELSE IF(QSIE.GE.0.D0)THEN
         IF(NRES.GT.0) WRITE(NRES,108)
108      FORMAT(/2X,' ZGOUBI  type  entrance  fringe  field'/)
C108      FORMAT(/2X,' Champ DE FUITE ENTREE ZGOUBI',/)
         IQSIE=1
      ENDIF
C
      UMEGA = A(NOEL,25)
      THETA = A(NOEL,26)
      R1    = A(NOEL,27)
      U1    = A(NOEL,28)
      U2    = A(NOEL,29)
      R2    = A(NOEL,30)
C
C     ** POUR LES BESOINS DE LA SIMPLE PRECISION :
      IF(R1*R1 .GE. 1.D10) U1 = -1.D6
      IF(R2*R2 .GE. 1.D10) U2 =  1.D6

      IF(NRES.GT.0) THEN
        WRITE(NRES,106) LAMBDE
106     FORMAT (/,5X,'Entrance  EFB',/,10X,
     >  'Fringe  field :  LAMBDA =',F7.2,' CM')
C106     FORMAT (/,5X,'FACE  D''ENTREE',/,10X,
C     >  'Champ FUITE : LAMBDA =',F7.2,' CM')
        IF(LAMBDE .EQ. 0.D0) THEN
          WRITE(NRES,*) '            ( No  fringe  field ) '
C          WRITE(NRES,*) '            ( PAS DE Champ DE FUITE ) '
        ELSE
          IF(IQSIE .EQ. 0) THEN
            WRITE(NRES,127) (CI(I),I=1,6),NCOEFE,SHIFTE
127         FORMAT (10X,' COEFFICIENTS :',6F10.5
     1      ,/,10X,' COEFFICIENTS 1 A ',I1,' PRIS EN COMPTE'
     2      ,/,10X,' DECALAGE  FACE  MAGNETIQUE = ',G12.4,' CM',/)
          ELSEIF(IQSIE .EQ. 1) THEN
            WRITE(NRES,166) QSIE
166         FORMAT(10X,'QSI ENTREE =',F7.2,' CM')
          ENDIF
        ENDIF
        WRITE(NRES,103) UMEGA,THETA,R1,U1,U2,R2
        IF(R1*R2 .EQ. 0.D0) WRITE(NRES,123)
 123    FORMAT(10('*'),' ATTENTION !    R1 or R2 is zero',10('*'))
      ENDIF
C
      UMEG = UMEGA * RAD
      TETA = THETA * RAD
      UT = UMEG - TETA
      SINO = SIN(UT)
      COSO = COS(UT)
      TANO = SINO / COSO
C
C  CALCUL DES PARAMETRES DE FACE  :
C  AXE X DU REFERENTIEL  = PARALLELE  DIRCTN ACENT, SENS DES RAYONS CROISSANTS
C  AXE Y DU REFERENTIEL  = ORTHOGONAL DIRCTN ACENT, DIRIGE VERS FACE D'ENTREE
C  ORIGINE DU REFERENTIEL= A L'INTERSECTION DE RM ET ACENT
C
C  PROJECTIONS DU RAYON MOYEN SUR LES AXES X,Y
      XB = RM * ( COS(UMEG) - 1D0 )
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
C  PENTE DE LA PARTIE LINEAIRE DU ChampDE FUITE 'ZGOUBI'
      IF(IQSIE .EQ. 1) THEN
         IF(QSIE +LAMBDE .GT.0.D0) THEN
            LDRTE=-1.D0/(QSIE +LAMBDE )
         ELSE
            LDRTE=99999.D0
         ENDIF
         PARAE  = LAMBDE**2 - QSIE**2
      ENDIF
C     WRITE(NRES,111) XA ,YA ,XB ,YB ,XC ,YC ,XD ,YD ,XE ,YE
C111  FORMAT(10G12.3)
C
C Champ DE FUITE SORTIE
      LAMBDS = A(NOEL,31)
      QSIS   = A(NOEL,32)
      NCOEFS = NINT(A(NOEL,33))
      DO 228 I=1,6
 228    CS(I) = A(NOEL,33+I)
      SHIFTS = A(NOEL,40)
      SHIFTS=0D0

      NCFS1=NCOEFS+1
      IF(QSIS.LT.0.D0) THEN
         IF(NRES.GT.0) WRITE(NRES,117)
117      FORMAT(/2X,' Exponential  exit  fringe  field',/)
         QSIS=0.D0
         IQSIS=0
      ELSE IF(QSIS.GE.0.D0)THEN
         IF(NRES.GT.0) WRITE(NRES,118)
118      FORMAT(/2X,' ZGOUBI  type  exit  fringe  field',/)
         IQSIS=1
      ENDIF
C
      UMEGAS = A(NOEL,41)
      THETAS = A(NOEL,42)
      R1S    = A(NOEL,43)
      U1S    = A(NOEL,44)
      U2S    = A(NOEL,45)
      R2S    = A(NOEL,46)
C
C     ** POUR LES BESOINS DE LA SIMPLE PRECISION :
      IF(R1S*R1S .GE. 1.D10) U1S = -1.D6
      IF(R2S*R2S .GE. 1.D10) U2S =  1.D6
C
      IF(NRES.GT.0) THEN
        WRITE(NRES,116) LAMBDS
116     FORMAT (/,5X,'Exit  EFB',/,10X,
     >  'Fringe  field :  Lambda =',F7.2,' CM')
        IF(LAMBDS .EQ. 0.D0) THEN
          WRITE(NRES,*) '            (No  Fringe  field ) '
        ELSE
          IF(IQSIS .EQ. 0) THEN
            WRITE(NRES,127) (CS(I),I=1,6),NCOEFS,SHIFTS
          ELSEIF(IQSIS .EQ. 1) THEN
            WRITE(NRES,167) QSIS
167         FORMAT(10X,'QSI  exit =',F7.2,' CM')
          ENDIF
        ENDIF
        WRITE(NRES,103) UMEGAS,THETAS,R1S,U1S,U2S,R2S
        IF(R1S*R2S .EQ. 0.D0) WRITE(NRES,123)
      ENDIF
C
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
      IF(IQSIS .EQ. 1) THEN
         IF(QSIS+LAMBDS.GT.0.D0) THEN
            LDRTS=-1.D0/(QSIS+LAMBDS)
         ELSE
            LDRTS=99999.D0
         ENDIF
         PARAS = LAMBDS**2 - QSIS**2
      ENDIF
C     WRITE(NRES,111) XAS,YAS,XBS,YBS,XCS,YCS,XDS,YDS,XES,YES
C
C Champ DE FUITE FACE EXTERNE (SEULEMENT SI OPTION 3 FACES)
      IF(NBFACE .EQ. 3) THEN
       LAMBD3 = A(NOEL,47)
       QSI3   = A(NOEL,48)
       NCOEF3 = NINT(A(NOEL,49))
       DO 229 I=1,6
 229     C3(I) = A(NOEL,49+I)
       SHIFT3 = A(NOEL,56)
       SHIFT3=0D0

       NCF31=NCOEF3+1
       IF(QSI3.LT.0.D0) THEN
          IF(NRES.GT.0) WRITE(NRES,137)
137       FORMAT(/2X,' Exponential  lateral  fringe  field',/)
          QSI3 =0D0
          IQSI3=0
       ELSE IF(QSI3.GE.0.D0)THEN
          IF(NRES.GT.0) WRITE(NRES,138)
138       FORMAT(/2X,' ZGOUBI  type  lateral  fringe  field',/)
          IQSI3=1
       ENDIF
C
       UMEGA3 = A(NOEL,57)
       THETA3 = A(NOEL,58)
       R13    = A(NOEL,59)
       U13    = A(NOEL,60)
       U23    = A(NOEL,61)
       R23    = A(NOEL,62)
       RM3    = A(NOEL,63)

       IF(NRES.GT.0) THEN
         WRITE(NRES,136) LAMBD3
136      FORMAT (/,5X,'Lateral  EFB',/,10X,
     >   'Fringe  field  : LAMBDA =',F7.2,' CM')
         IF(LAMBD3 .EQ. 0.D0) THEN
           WRITE(NRES,*) '            (No  fringe  field) '
         ELSE
           IF(IQSI3 .EQ. 0) THEN
             WRITE(NRES,127) (C3(I),I=1,6),NCOEF3,SHIFT3
           ELSEIF(IQSI3 .EQ. 1) THEN
             WRITE(NRES,135) QSI3
135          FORMAT(10X,'QSI lateral  =',F7.2,' CM')
           ENDIF
         ENDIF
         WRITE(NRES,113) RM3
113      FORMAT(20X,' Face centered on direction ACENT+OMEGA, A ',
     >   G12.4,' CM'/)
         WRITE(NRES,103) UMEGA3,THETA3,R13,U13,U23,R23
         IF(R13*R23 .EQ. 0.D0) WRITE(NRES,123)
       ENDIF
C
       UMEG3 = UMEGA3 * RAD
       TETA3 = THETA3 * RAD
       UT3= UMEG3- TETA3
       IF(ABS(UT3).GE. .5D0*PI) THEN
         IF(NRES.GT.0) WRITE(NRES,139)
139      FORMAT(/,20X,' FACE 3, ATTENTION : ABS(UT3) > 90 DEG. =>',
     >   ' Value is ambiguous for calculation of field map')
       ENDIF
C
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
       IF(IQSI3 .EQ. 1) THEN
          IF(QSI3+LAMBD3.GT.0.D0) THEN
             LDRT3=-1.D0/(QSI3+LAMBD3)
          ELSE
             LDRT3=99999.D0
          ENDIF
          PARA3 = LAMBD3**2 - QSI3**2
       ENDIF
C      WRITE(NRES,111) XA3,YA3,XB3,YB3,XC3,YC3,XD3,YD3,XE3,YE3
      ENDIF
C     ... NFACE=3
C
C CALCUL LA CARTE DE Champ, POUR IXMA ANGLES ET JYMA RAYONS
C
C     ****DEBUT DE BOUCLES SUR ANGLES ET RAYONS
      DO  1  J = 1,IXMA
      NN = NN + 1
      XH(NN) = (NN - 1) * ATS * RAD
      ZETA = ACN - XH(NN)
      DO  2  I = 1,JYMA
C
C     ****COORDONNEES DU POINT COURANT
      X = YH(I) * COS(ZETA) - RM
      Y = YH(I) * SIN(ZETA)
C     ****AX (CX) = COTE X DE LA PROJECTION DE A (C) SUR LA DROITE (MX)
C     (M=POINT COURANT), ORTHOGONALEMENT A LA DROITE (ABC)
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
      LIEU=-1
      IFACE = 10
      IF     ( X.LT.CX .AND. X.LT.AX )  THEN
C        REGION DE COURBURE R1
         RR = R1
         XO = XD
         YO = YD
         R = ABS(RR)
         DPE = R - SQRT((X - XO)**2 + (Y - YO)**2 )
         IF(DPE*RR .LT. 0.D0) IFACE = 1
      ELSE IF( X.GE.CX .AND. X.LE.AX )  THEN
C        REGION LINEAIRE
         XO = SICO * (Y - YB) + XB * SIN2 + X * COS2
         YO = SICO * (X - XB) + YB * COS2 + Y * SIN2
         YL = YB + ((X - XB) * TANO )
C        ****YL = COTE Y DE LA PROJECTION DE B SUR LA DROITE (MY)
C        (M=POINT COURANT), PARALLELEMENT A LA DROITE (ABC)
         DPE = SQRT((X - XO)**2 + (Y - YO)**2 )
         IF( Y .LE. YL .OR. DPE .LE. 1.D-6 ) THEN
C           COTE INTERIEUR / FACE ENTREE
            LIEU=4
         ELSE
C           COTE EXTERIEUR / FACE ENTREE
            LIEU=14
            IFACE = 1
         ENDIF
      ELSE IF( X.GT.CX .AND. X.GT.AX )  THEN
C        REGION DE COURBURE R2
         RR = R2
         XO = XE
         YO = YE
         R = ABS(RR)
         DPE = R - SQRT((X - XO)**2 + (Y - YO)**2 )
         IF(DPE*RR .LT. 0.D0) IFACE = 1
      ELSE
C        ERREUR  DE  DONNEES  FACE  ENTREE
         IF(NRES.GT.0) WRITE(NRES,104)
104      FORMAT(/,5X,10('*'),' ERROR  on  parameters  entrance  face',/)
         GOTO  15
      ENDIF
C
C POSITION DU POINT (X,Y) / FACE SORTIE
         LIEUS=-1
         IFACES = 20
         IF     ( X.LT.CXS .AND. X.LT.AXS )  THEN
C           REGION DE COURBURE R1
            RRS = R1S
            XO = XDS
            YO = YDS
            R = ABS(RRS)
            DPS = R - SQRT((X - XO)**2 + (Y - YO)**2 )
            IF(DPS*RRS .GT. 0.D0) IFACES = 2
         ELSE IF( X.GE.CXS .AND. X.LE.AXS )  THEN
C           REGION LINEAIRE
            XO = SICOS * (Y - YBS) + XBS * SIN2S + X * COS2S
            YO = SICOS * (X - XBS) + YBS * COS2S + Y * SIN2S
            YL = YBS + ((X - XBS) * TANOS )
            DPS = SQRT((X - XO)**2 + (Y - YO)**2 )
            IF( Y .GE. YL .OR. DPS .LE.1.D-6 ) THEN
C              COTE INTERIEUR / FACE SORTIE
               LIEUS=4
            ELSE
C              COTE EXTERIEUR / FACE SORTIE
               LIEUS=14
               IFACES = 2
            ENDIF
         ELSE IF( X.GT.CXS .AND. X.GT.AXS )  THEN
C           REGION DE COURBURE R2
            RRS = R2S
            XO = XES
            YO = YES
            R = ABS(RRS)
            DPS = R - SQRT((X - XO)**2 + (Y - YO)**2 )
            IF(DPS*RRS .GT. 0.D0) IFACES = 2
         ELSE
            IF(NRES.GT.0) WRITE(NRES,114)
114         FORMAT(/,5X,10('*'),' ERROR  on  parameters  exit  face',/)
            GOTO  15
         ENDIF
C
C POSITION DU POINT (X,Y) / FACE EXTERNE (SEULEMENT SI OPTION 3 FACES)
      IF(NBFACE .EQ. 3) THEN
          LIEU3=-1
          IFACE3 = 30
          IF     ( X.LT.CX3 .AND. X.LT.AX3 )  THEN
C           REGION DE COURBURE R1
            RR3 = R13
            XO = XD3
            YO = YD3
            R = ABS(RR3)
            DP3 = R - SQRT((X - XO)**2 + (Y - YO)**2 )
            IF(DP3*RR3 .GT. 0.D0) IFACE3 = 3
          ELSE IF( X.GE.CX3 .AND. X.LT.AX3 )  THEN
C           REGION LINEAIRE
            XO = SICO3 * (Y - YB3) + XB3 * SIN23 + X * COS23
            YO = SICO3 * (X - XB3) + YB3 * COS23 + Y * SIN23
            YL = YB3 + ((X - XB3) * TANO3 )
            DP3 = SQRT((X - XO)**2 + (Y - YO)**2 )
            IF( Y .GE. YL .OR. DP3 .LE. 1.D-6 ) THEN
C              COTE INTERIEUR / FACE EXTERNE
               LIEU3=4
            ELSE
C              COTE EXTERIEUR / FACE EXTERNE
               LIEU3=14
               IFACE3 = 3
            ENDIF
          ELSE IF( X.GE.CX3 .AND. X.GE.AX3 )  THEN
C           REGION DE COURBURE R2
            RR3 = R23
            XO = XE3
            YO = YE3
            R = ABS(RR3)
            DP3 = R - SQRT((X - XO)**2 + (Y - YO)**2 )
            IF(DP3*RR3 .GT. 0.D0) IFACE3 = 3
          ELSE
            IF(NRES.GT.0) THEN
              WRITE(NRES,124)
124           FORMAT(/,5X,'*** ERROR  on  parameters lateral  face',/)
              WRITE(NRES,134) X,AX3,CX3,XA3,YA3,XC3,YC3
134           FORMAT(/,5X,7G12.4)
            ENDIF
            GOTO  15
          ENDIF
        ENDIF
C
CE 'POINT DE MESURE' S'UTILISE AVEC UNE TRAJECTOIRE. IL PERMET DE CONNAITRE
C  LA DISTANCE DU POINT COURANT X,Y A CHACUNE DES 3 FACES, RESPECTIVEMENT
C  DPE(ENTREE), DPS(SORTIE), DP3(FACE 'EXTERNE')
C      IF(I .EQ. JYMA/2)
C     X WRITE(NRES,8002) IFACE,IFACES,IFACE3
C 8002 FORMAT(/,' IFACE,S,3 :', 3I3)
C     X WRITE(NRES,8002) LIEU,DPE,RR,LIEUS,DPS,RRS,LIEU3,DP3,RR3
C 8002 FORMAT(/,' LIEU,DPE,RR  E/S/3 : ',/,3(I3,2G12.4))
C
C
C        ** CALCULE LES COEFFS. DE Champ DE FUITE ASSOCIES AUX FACES

         LO=ABS(DPE)
         IF(IFACE .EQ. 1) THEN
C           ****EXTERIEUR AIMANT / FACE D'ENTREE
            IF(LO.LT.LAMBDE) THEN
               IF(IQSIE .EQ. 0) THEN
C                 ****Champ DE FUITE EXPONENTIEL
                  IF(LIEU .EQ. -1) THEN
C                    ****ZONE R1 OU R2
                     DPNORM=(-DPE*RR/ABS(RR)+SHIFTE)/LAMBDE
                  ELSE
C                    ****ZONE LINEAIRE U1 OU U2
                     DPNORM=( DPE           +SHIFTE)/LAMBDE
                  ENDIF
                  PH=CI(1)
                  DO 18 ICOEF=2,NCFE1
18                   PH=PH+CI(ICOEF)*DPNORM**(ICOEF-1)
                  FH=1.D0/(1.D0+EXP(PH))
               ELSEIF(IQSIE .EQ. 1) THEN
C                 ****Champ DE FUITE 'ZGOUBI'
                  IF     ( LO.LE.QSIE )  THEN
                     FH =    (LO * LDRTE) + 0.5D0
                  ELSE IF( LO.GT.QSIE )  THEN
                     FH =    ((LO - LAMBDE)**2 ) * 0.5D0 / PARAE
                  ENDIF
               ENDIF
            ELSEIF(LO.GE.LAMBDE) THEN
               IF(LO .EQ. 0.D0) THEN
                 FH=1.D0
               ELSE
                 FH=0.D0
               ENDIF
            ENDIF
         ELSEIF(IFACE .EQ. 10) THEN
C           ****INTERIEUR AIMANT / FACE D'ENTREE
            IF(LO.LT.LAMBDE) THEN
               IF(IQSIE .EQ. 0.D0) THEN
C                 ****Champ DE FUITE EXPONENTIEL
                  IF(LIEU .EQ. -1) THEN
C                    ****ZONE R1 OU R2
                     DPNORM=(-DPE*RR/ABS(RR)+SHIFTE)/LAMBDE
                  ELSE
C                    ****ZONE LINEAIRE U1 OU U2
                     DPNORM=(-DPE           +SHIFTE)/LAMBDE
                  ENDIF
                  PH=CI(1)
                  DO 19 ICOEF=2,NCFE1
19                   PH=PH+CI(ICOEF)*DPNORM**(ICOEF-1)
                  FH=1.D0/(1.D0+EXP(PH))
               ELSEIF(IQSIE .EQ. 1) THEN
C                 ****Champ DE FUITE 'ZGOUBI'
                  IF     ( LO.LE.QSIE )  THEN
                     FH = .5D0-(LO * LDRTE)
                  ELSE IF( LO.GT.QSIE )  THEN
                     FH = 1.D0-((LO - LAMBDE)**2 ) * 0.5D0 / PARAE
                  ENDIF
               ENDIF
            ELSEIF(LO.GE.LAMBDE) THEN
               FH=1.D0
            ENDIF
         ENDIF
         FHE = FH

         LO=ABS(DPS)
         IF(IFACES .EQ. 2) THEN
C           ****EXTERIEUR AIMANT / FACE DE SORTIE
            IF(LO.LT.LAMBDS) THEN
               IF(IQSIS .EQ. 0) THEN
C                 ****Champ DE FUITE EXPONENTIEL
C                 SIGNE= RRS/ABS(RRS)
C                 DPNORM=(DPS*SIGNE+SHIFTS)/LAMBDS
                  IF(LIEUS .EQ. -1) THEN
C                    ****ZONE R1 OU R2
                     DPNORM=( DPS*RRS/ABS(RRS)+SHIFTS)/LAMBDS
                  ELSE
C                    ****ZONE LINEAIRE U1 OU U2
                     DPNORM=( DPS           +SHIFTS)/LAMBDS
                  ENDIF
                  PH=CS(1)
                  DO 33 ICOEF=2,NCFS1
33                   PH=PH+CS(ICOEF)*DPNORM**(ICOEF-1)
                  FH=1.D0/(1.D0+EXP(PH))
               ELSEIF(IQSIS .EQ. 1) THEN
C                 ****Champ DE FUITE 'ZGOUBI'
                  IF     ( LO.LE.QSIS )  THEN
                     FH =    (LO * LDRTS) + 0.5D0
                  ELSE IF( LO.GT.QSIS )  THEN
                     FH =    ((LO - LAMBDS)**2 ) * 0.5D0 / PARAS
                  ENDIF
               ENDIF
            ELSEIF(LO.GE.LAMBDS) THEN
               IF(LO .EQ. 0.D0) THEN
                 FH=1.D0
               ELSE
                 FH=0.D0
               ENDIF
            ENDIF
         ELSE IF(IFACES .EQ. 20) THEN
C           ****INTERIEUR AIMANT / FACE SORTIE
            IF(LO.LT.LAMBDS) THEN
               IF(IQSIS .EQ. 0) THEN
C                 ****Champ DE FUITE EXPONENTIEL
C                 SIGNE= RRS/ABS(RRS)
C                 DPNORM=(DPS*SIGNE+SHIFTS)/LAMBDS
                  IF(LIEUS .EQ. -1) THEN
C                    ****ZONE R1 OU R2
                     DPNORM=( DPS*RRS/ABS(RRS)+SHIFTS)/LAMBDS
                  ELSE
C                    ****ZONE LINEAIRE U1 OU U2
                     DPNORM=(-DPS           +SHIFTS)/LAMBDS
                  ENDIF
                  PH=CS(1)
                  DO 27 ICOEF=2,NCFS1
27                   PH=PH+CS(ICOEF)*DPNORM**(ICOEF-1)
                  FH=1.D0/(1.D0+EXP(PH))
               ELSEIF(IQSIS .EQ. 1) THEN
C                 ****Champ DE FUITE 'ZGOUBI'
                  IF     ( LO.LE.QSIS )  THEN
                     FH = .5D0-(LO * LDRTS)
                  ELSE IF( LO.GT.QSIS )  THEN
                     FH = 1.D0-((LO - LAMBDS)**2 ) * 0.5D0 / PARAS
                  ENDIF
               ENDIF
            ELSEIF(LO.GE.LAMBDS) THEN
               FH=1D0
            ENDIF
         ENDIF
         FHS = FH

         IF(NBFACE .EQ. 3) THEN
         LO=ABS(DP3)
         IF(IFACE3 .EQ. 3) THEN
C           ****EXTERIEUR AIMANT / FACE EXTERNE
            IF(LO.LT.LAMBD3) THEN
               IF(IQSI3 .EQ. 0) THEN
C                 ****Champ DE FUITE EXPONENTIEL
C                 SIGNE= RR3/ABS(RR3)
C                 DPNORM=(DP3*SIGNE+SHIFT3)/LAMBD3
                  IF(LIEU3 .EQ. -1) THEN
C                    ****ZONE R1 OU R2
                     DPNORM=( DP3*RR3/ABS(RR3)+SHIFT3)/LAMBD3
                  ELSE
C                    ****ZONE LINEAIRE U1 OU U2
                     DPNORM=( DP3           +SHIFT3)/LAMBD3
                  ENDIF
                  PH=C3(1)
                  DO 29 ICOEF=2,NCF31
29                   PH=PH+C3(ICOEF)*DPNORM**(ICOEF-1)
                  FH=1.D0/(1.D0+EXP(PH))
               ELSEIF(IQSI3 .EQ. 1) THEN
C                 ****Champ DE FUITE 'ZGOUBI'
                  IF     ( LO.LE.QSI3 )  THEN
                     FH =    (LO * LDRT3) + 0.5D0
                  ELSE IF( LO.GT.QSI3 )  THEN
                     FH =    ((LO - LAMBD3)**2 ) * 0.5D0 / PARA3
                  ENDIF
               ENDIF
            ELSEIF(LO.GE.LAMBD3) THEN
               IF(LO .EQ. 0.D0) THEN
                 FH=1D0
               ELSE
                 FH=0D0
               ENDIF
            ENDIF
         ELSE IF(IFACE3 .EQ. 30) THEN
C           ****INTERIEUR AIMANT / FACE EXTERNE
            IF(LO.LT.LAMBD3) THEN
               IF(IQSI3 .EQ. 0) THEN
C                 ****Champ DE FUITE EXPONENTIEL
C                 SIGNE= RR3/ABS(RR3)
C                 DPNORM=(DP3*SIGNE+SHIFT3)/LAMBD3
                  IF(LIEU3 .EQ. -1) THEN
C                    ****ZONE R1 OU R2
                     DPNORM=( DP3*RR3/ABS(RR3)+SHIFT3)/LAMBD3
                  ELSE
C                    ****ZONE LINEAIRE U1 OU U2
                     DPNORM=(-DP3           +SHIFT3)/LAMBD3
                  ENDIF
                  PH=C3(1)
                  DO 34 ICOEF=2,NCF31
34                   PH=PH+C3(ICOEF)*DPNORM**(ICOEF-1)
                  FH=1.D0/(1.D0+EXP(PH))
               ELSEIF(IQSI3 .EQ. 1) THEN
C                 ****Champ DE FUITE 'ZGOUBI'
                  IF     ( LO.LE.QSI3 )  THEN
                     FH = .5D0-(LO * LDRT3)
                  ELSE IF( LO.GT.QSI3 )  THEN
                     FH = 1.D0-((LO - LAMBD3)**2 ) * 0.5D0 / PARA3
                  ENDIF
               ENDIF
            ELSEIF(LO.GE.LAMBD3) THEN
               FH=1D0
            ENDIF
         ENDIF
         FH3 = FH
         ENDIF

      FH = FHE * FHS
      IF(NBFACE .EQ. 3) FH = FH * FH3

C     CE 'POINT DE MESURE' PERMET DE SUIVRE LE
C     PARCOURS A TRAVERS LES FACES, LE LONG D'UN RAYON.
C      IF(I .EQ. JYMA/2)
C     X WRITE(NRES,8001) FHE,FHS,FH3,FH,X,Y,J,I
C 8001 FORMAT(' .....FHE,FHS,FH3,FH  ,X,Y, ANG,RAY :',6G12.4,2I5)
C
      YHI= (YH(I)-RM)/RM
      HC(ID,NN,I,1,IMAP)=
     >      FH*HNORM*(1.D0+(COEFN+(COEFB+COEFG*YHI)*YHI)*YHI)
C
C     ****FIN DE BOUCLES SUR ANGLES ET RAYONS
    2 CONTINUE
    1 CONTINUE

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
109         FORMAT(/,20X,' Perturbation  to  field  map, linear'
     >      ,'EN ANGLE ET',/,20X,' CENTREE SUR THETA =',F7.2
     >      ,' DEG.,  AVEC DB/B0(TOTALE) = ',F10.6)
            J0=INT(ZETA/ATS)
            DZETA=ZETA-J0*ATS
            DBSB=DBSB /IXMA
            PH=1.D0-DBSB*DZETA/ATS
            DO 3 J=1,IXMA
              DO 3 I=1,JYMA
                HC(ID,J ,I ,1,IMAP)=HC(ID,J ,I, 1,IMAP)*(PH+DBSB*(J-J0))
3           CONTINUE
         ELSEIF(NBSHIM .EQ. -2) THEN
C           ****PERTURBATION LINEAIRE EN RAYON
            RR   = A(NOEL,ND+1)
            DBSB = A(NOEL,ND+2)
            IF(NRES.GT.0) WRITE(NRES,110) RR  ,DBSB
110         FORMAT(/,20X,' Perturbation  to  field  map, linear'
     >      ,'EN RAYON ET',/,20X,' CENTREE SUR R =',F7.2
     >      ,' CM,   AVEC DB/B0(TOTALE) = ',F10.6)
            J0=INT(RR  /RMS)
            DRR  =RR  -J0*RMS
            DBSB=DBSB /JYMA
            PH=1.D0-DBSB*DRR  /RMS
            DO 4 J=1,IXMA
              DO 4 I=1,JYMA
                HC(ID,J ,I, 1,IMAP)=HC(ID,J ,I, 1,IMAP)*(PH+DBSB*(I-J0))
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

  100 FORMAT(20X,'AIMANT  PRINCIPAL',///,
     1 11X,'ANGLES : A.TOTAL =',F6.2,' degrees',5X,'A.CENTRAL =',
     2 F6.2,' degrees',// ,11X,'RM =',F7.2,' cm',5X,'RMIN =',F7.2,' cm',
     3 5X,'RMAX =',F7.2,' cm',//,
     5 11X,'HNORM =',F8.4,' kGauss',5X,'COEF.N =',F9.5,5X,'COEF.B =',
     6 F9.5,5X,'COEF.G=',F9.5)
C  103 FORMAT(10X,7HOMEGA =,F7.2,5X,17HANGLE  DE  FACE =,F7.2,// ,
 103  FORMAT(10X,'OMEGA =',F7.2,' deg.',
     > 5X,'Wedge  angle  =',F7.2,' deg.',/,
     1 11X,'Radius 1 =',1P,G10.2,' cm,'// ,
     2 11X,'Straight  segment 1 =',G10.2,' cm',// ,
     3 11X,'Straight  segment 2 =',G10.2,' cm',// ,
     4 11X,'Radius 2 =',G10.2,' cm')

      BMIN = -1.1D8
      BMAX = 1.1D8
      BNORM = 1.2D8
      ABMI = -1.3D8
      RBMI = -1.4D8
      ZBMI = -1.5D8
      ABMA = 1.3D8
      RBMA = 1.4D8
      ZBMA = 1.5D8

      RETURN
      END
