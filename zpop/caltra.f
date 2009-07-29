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
      SUBROUTINE CALTRA(N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,20)
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/REBELO/ NRBLT,IPASS,KWRI,NNDES,STDVM
      CHARACTER FAM*8,LBF*8,KLEY*10,LABEL*8
      INCLUDE 'MXFS.H'
      COMMON/SCALT/ FAM(MXF),LBF(MXF,2),KLEY,LABEL(MXL,2)
      CHARACTER*80 TITRE
      COMMON/TITR/ TITRE 

      PARAMETER(MPOL=10)
 
      DIMENSION ND(MXL)

      LOGICAL TOMANY

      LOGICAL READAT

      COMMON/VEN/Y1,S1,D,DQ,DS,DM,DSS,DPU,TETA1,
     >           SMIN,SMAX,YMIN,YMAX,SB0,YB0
      INTEGER BLACK, BLUE, CYAN

      CHARACTER*32 TXT

      DATA BLACK, BLUE, CYAN / 0,0,0 /

      PARAMETER (ZERO = 0.D0)

C This INCLUDE must stay located right before the first statement
      INCLUDE 'LSTKEY.H'

      PI = 4.D0*ATAN(1.D0) 
      DEG2RD = PI / 180.D0
      CM2M = 1.D-2

      DEV2 = 0.D0

      READAT = (N .EQ. 0)
     
      IF(READAT) CALL PRDATA(LABEL,NB)

      IF(READAT) THEN
        READ(NDAT,FMT='(A)') TITRE
        WRITE(6,FMT='(/,A,/)') TITRE
      ENDIF

      TOMANY = .FALSE.
      NOEL = 0
      S=0.D0
      Y=0.D0
      Z=0.D0
      TETA1=0.D0
      TETAZ=0.D0
      TETREF = 0.D0

CC(1X,A,1P,2G12.4,1X,2G12.4,1X,3G12.4)')
      IF(N.EQ.1) WRITE(6,FMT='(/,T2,A,T9,A,T21,A,T34,A,T46,A,T58,A,
     > T70,A,T82,A)')  'NAME','S_in','Y_in',
     >  'S_out','Y_out','Theta_in','Theta_out','Width :'

      GOTO 990

 998  CONTINUE

 990  CONTINUE

      IF(READAT) THEN
        READ(NDAT,*) KLEY
          write(*,*) ' zpop_synop ',kley
        DO 188 IKLE=1,MXKLE
          IF(KLEY .EQ. KLE(IKLE)) THEN
            NOEL = NOEL+1
            IF( NOEL .EQ. MXL+1) THEN
              TOMANY=.TRUE.
              WRITE(6,*) ' QUEUE DATA WERE SKIPPED: too many elements'
              WRITE(6,*) ' (number of elements should not exceed '
     >                                                       ,MXL,').'
              GOTO 99
            ENDIF
            IQ(NOEL) =  IKLE
            GOTO 187
          ENDIF
 188    CONTINUE

      ELSE

        NOEL = NOEL+1
        IKLE = IQ(NOEL)
        KLEY = KLE(IKLE)

      ENDIF
 
 187  CONTINUE

      COEF=COS(TETAZ)
      SP=S
      YP=Y 

C Go to "GOTO(...) IKLE" ---------------------------
      GOTO 1001
  999 CONTINUE
C---------------------------------------------------

      CALL FBGTXT
      WRITE(6,*) '  ******* WARNING : unknown key ',KLEY
      GOTO 998

 995  CONTINUE
      CALL FBGTXT
      IF(N.EQ.1) WRITE(6,FMT='(/,T2,A,T9,A,T21,A,T34,A,T46,A,T58,A,
     > T70,A,T82,A)')  'NAME','S_in','Y_in',
     >  'S_out','Y_out','Theta_in','Theta_out','Width'
      WRITE(6,200) KLEY
 200  FORMAT(/,10X,' Data list :  stopped upon key  ',A8,//,90(1H*))
      RETURN
 
C-----  DRIFT, ESL. ESPACE LIBRE
 1    CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1)
      IF(READAT) GOTO 998
      AL = A(NOEL,1) *COEF * CM2M 
      IF(AL.LT.0.D0) THEN
        WRITE(6,FMT='(/,A,/)')
     >    ' ** WARNING : cannot handle negative drift'
        GOTO 998
      ENDIF
      XI=0.D0
      PHI=0.D0
      CALL STORIG(NOEL,S1,Y1,Z,TETA1,XI,PHI)
C      CALL SYNAXE(N,TETA1,S1,Y1,TETA,S,Y,AL,DS,BLACK,BLACK)
      REDUC=AL/2.D0
      WIDTH = 0.D0   !!!!!!!DSS * REDUC  
      CALL SYNBOX(N,TETA1,S1,Y1,TETA,S,Y,AL,WIDTH)
      TXT = 'DRIF'
      IF(N.EQ.1) THEN
        WRITE(6,FMT='(1X,A,1P,2G12.4,1X,2G12.4,1X,3G12.4)')
     >           TXT,SP,YP,S,Y,TETA,TETA1,WIDTH
      ENDIF
      GOTO 998
C----- AIMANT. DIPOLE AVEC CARTE DE CHAMP POLAIRE CALCULEE (PGM CARLA)
 2    CONTINUE
      IF(READAT) CALL RAIMAN(NDAT,NOEL,MXL,A,ND(NOEL))
      IF(READAT) GOTO 998
      GOTO 998
C----- QUADRUPO
 3    CONTINUE
      IF(READAT) CALL RQSOD(NDAT,NOEL,MXL,A,ND(NOEL))
      IF(READAT) GOTO 998
 301  CONTINUE
      AL=A(NOEL,10)*COEF  * CM2M 
C      AK=A(NOEL,12)/A(NOEL,11)/BORO *1.D4
C      shift of origin (due to upstream fringe field extent) : 
      REDUC=AL/2.D0
      WIDTH = DQ*REDUC
      TXT='QUAD'
 310  CONTINUE

      XI = -A(NOEL,20)  * CM2M  
      IF(A(NOEL,21).LE.0.D0) XI = 0.D0

      PHI=0.D0
C       call fbgtxt
C         write(*,*) ' quadru  S,Y,Z,TETA1,XI ',S,Y,Z,TETA1,XI
      CALL STORIG(NOEL,S1,Y1,Z,TETA1,XI,PHI)
      CALL SYNBOX(N,TETA1,S1,Y1,TETA,S,Y,AL,WIDTH)
      IF(N.EQ.1) THEN
        WRITE(6,FMT='(1X,A,1P,2G12.4,1X,2G12.4,1X,3G12.4)')
     >           TXT,SP,YP,S,Y,TETA,TETA1,WIDTH
      ELSEIF(N.EQ.2) THEN
C        CALL TRTXT(S,Y,LABEL(NOEL,1),4,1)
C        CALL TRTXT(S,Y,TXT,4,1)
      ENDIF
      GOTO 998
C----- B SEXTUPOLAIRE  ET DERIVEES CALCULES EN TOUT POINT (X,Y,Z)
4     CONTINUE
      IF(READAT) CALL RQSOD(NDAT,NOEL,MXL,A,ND(NOEL))
      IF(READAT) GOTO 998
      GOTO 301
      GOTO 998
C----- RECHERCHE DU PLAN IMAGE ET DIMENSIONS D'IMAGE HORIZONTAUX
C----- FAISCEAU MONOCHROMATIQUE
 5    CONTINUE
      GOTO 998
C----- RECHERCHE DU PLAN IMAGE ET DIMENSIONS D'IMAGE HORIZONTAUX
C----- FAISCEAU POLYCHROME
 6    CONTINUE
      GOTO 998
C----- FAISCEAU - LISTE DU FAISCEAU AU POINT COURANT (DANS ZGOUBI.RES)
7     CONTINUE
      GOTO 998
C----- FAISCNL - LISTE DU FAISCEAU AU POINT COURANT (DANS ZGOUBI.FAI)
 8    CONTINUE
      IF(READAT)  READ(NDAT,*) TA(NOEL,1)
      IF(READAT) GOTO 998
      GOTO 998
C----- TRAVERSEE D'UNE CIBLE
9     CONTINUE
      IF(READAT) CALL RCIBLE
      IF(READAT) GOTO 998
C      CALL CIBLE
      GOTO 998
C----- DIMENSIONS DU FAISCEAU @ XI
10    CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1)
      IF(READAT) GOTO 998
      GOTO 998
C----- REBELOTE - RE-PASSE IPASS FOIS DANS LA STRUCTURE
11    CONTINUE
      IF(READAT) READ(NDAT,*) (A(NOEL,I),I=1,3)
      IF(READAT) GOTO 998
      GOTO 998
C----- QUADISEX. CHAMP CRENEAU B = B0(1+N.Y+B.Y2+G.Y3) PLAN MEDIAN
 12   CONTINUE
      KALC=1
      KUASEX=3
      IF(READAT) CALL RSIMB(ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- CHANGREF. TRANSLATION X,Y ET ROTATION DU REFERENTIEL COURANT
 13   CONTINUE
      IF(READAT) READ(NDAT,*) (A(NOEL,II),II=1,3)
      IF(READAT) GOTO 998
      XCE=A(NOEL,1) * CM2M
      YCE=A(NOEL,2) * CM2M
      ALE=A(NOEL,3) * DEG2RD
      XI = 0.D0
      PHI2 = -ALE 
      CALL STORIG(NOEL,S,Y,Z,TETA1-PHI2,XI,PHI2)
      AL = COEF * SQRT(XCE*XCE + YCE*YCE) 
      TETA=TETA1 + ATAN2(YCE,XCE)
      S=S1+AL*COS(TETA)
      Y=Y1+AL*SIN(TETA)
      CALL CALTRI(N,S,Y)
      TETA1 = TETA1 - PHI2
      GOTO 998
C----- SEXQUAD. CHAMP CRENEAU B = B0(N.Y+B.Y2+G.Y3) PLAN MEDIAN
14    CONTINUE
      IF(READAT) CALL RSIMB(ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- AIMANT TOROIDAL
15    CONTINUE
      IF(READAT) CALL RSIMB(ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- OBJET MONTE-CARLO = MU OU PI DE DECROISSANCE DU ETA
C          ( D'APRES B. MAYER , FEVRIER 1990 )
 17   CONTINUE
      IF(READAT) CALL ROBJTA
      IF(READAT) GOTO 998
      BORO = A(NOEL,1)
      GOTO 998
C----- COEFFICIENTS D'ABERRATION A L'ABSCISSE COURANTE
 18   CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1),A(NOEL,2)
      IF(READAT) GOTO 998
      GOTO 998
C----- COMPTE LES TRAJECTOIRES HORS LIMITES DE LA CHAMBRE
 19   CONTINUE
      IF(READAT) CALL RCOLLI
      IF(READAT) GOTO 998
      GOTO 998
C----- CHAMP Q-POLE SPECIAL PROJET SPES2
 20   CONTINUE
      KALC=1
      KUASEX=1
C ............... ADD   IF(READAT) CALL ............
C ...............RETURN
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- CARTEMES. CARTE DE CHAMP CARTESIENNE MESUREE DU SPES2
 21   CONTINUE
      KALC=2
      KUASEX=1
      IF(READAT) CALL RCARTE(1,2,ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- CHANGE LE FAISCEAU Y,T,Z,P EN -Y,-T,-Z,-P
C       ( EQUIVALENT A CHANGEMENT DE SIGNE DU CHAMP DIPOLAIRE }
 22   CONTINUE
C      CALL YMOINY
      GOTO 998
C----- B OCTUPOLAIRE ET DERIVEES CALCULES EN TOUT POINT (X,Y,Z)
 23   CONTINUE
      KALC=3
      KUASEX=4
      IF(READAT) CALL RQSOD(NDAT,NOEL,MXL,A,ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- OBJET
24    CONTINUE
      IF(READAT) CALL ROBJET
      BORO = A(NOEL,1)
      IF(READAT) GOTO 998
      GOTO 998
C----- OBJET DEFINI PAR MONTE-CARLO
 25   CONTINUE
      IF(READAT) CALL RMCOBJ
      IF(READAT) GOTO 998
      BORO = A(NOEL,1)
      GOTO 998
C----- DESINTEGRATION EN COURS DE VOL
26    CONTINUE
      IF(READAT) CALL RMCDES
      IF(READAT) GOTO 998
      GOTO 998
C----- HISTOGRAMME DE Y T Z P D
27    CONTINUE
      IF(READAT) CALL RHIST(NDAT,NOEL,MXL,A,TA)
      IF(READAT) GOTO 998
      GOTO 998
C----- TRANSMAT. TRANSFERT MATRICIEL AU SECOND ORDRE
 28   CONTINUE
      IF(READAT) CALL RTRANS
      IF(READAT) GOTO 998
C      CALL TRANSM
      GOTO 998
C----- AIMANT VENUS (CHAMP CONSTANT DANS UN RECTANBLE)
 29   CONTINUE
      KALC=1
      KUASEX=5
      IF(READAT) CALL RSIMB(ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- AIMANT PS170 (CHAMP CONSTANT DANS UN CERCLE)
 30   CONTINUE
      KALC=1
      KUASEX=6
      IF(READAT) CALL RSIMB(ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- TOSCA. CARTE DE CHAMP CARTESIENNE 2-D OU 3-D FABRIQUEE PAR TOSCA
 31   CONTINUE
      KALC=2
      IF(READAT) CALL RCARTE(KART,3,ND(NOEL))
      IF(READAT) GOTO 998
      IF(A(NOEL,22) .EQ. 1) THEN
        KUASEX=2
      ELSE
        KUASEX=7
      ENDIF
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- AUTOREF. CHGNT RFRNTIEL AUTMTIQ => PLAN NORMAL A LA TRAJ. DE REFERENCE
 32   CONTINUE
      IF(READAT) CALL RAUTOR
      IF(READAT) GOTO 998
C      CALL AUTORF
      GOTO 998
C----- COLLIMATOR ( PARTICULES OUT => IEX = -3 )
 33   CONTINUE
      IF(READAT) CALL RCOLLI
      IF(READAT) GOTO 998
C      CALL COLLIM
      GOTO 998
C----- MULTIPOLE
 34   CONTINUE
      IF(READAT) CALL RMULTI(NDAT,NOEL,MXL,A,MPOL,ND(NOEL))
      IF(READAT) GOTO 998
      XL=A(NOEL,2)
      AL = XL*COEF  * CM2M 
      REDUC=AL/2.D0
      WIDTH=DQ*  REDUC

C      shift of origin (due to upstream fringe field extent) : 
      XI = -A(NOEL,14)  * CM2M  
      IF(A(NOEL,15).LE.0.D0) XI = 0.D0

      B1  =A(NOEL,4)
      KPOS = A(NOEL,3+ND(NOEL))
      TXT = 'MULT'

C            call fbgtxt
C            write(*,*) ' caltra mult ND ',ND(NOEL),kpos

        IF    (KPOS.EQ.3) THEN
          XCE=A(NOEL,4+ND(NOEL)) * CM2M
          YCE=A(NOEL,5+ND(NOEL)) * CM2M
          ALE=A(NOEL,6+ND(NOEL)) 
          XIRF = 0.D0
          PHI2RF = -ALE 
          CALL STORIG(NOEL,S,Y,Z,TETA1-PHI2RF,XI,PHI2RF)
          ALSHFT = COEF * SQRT(XCE*XCE + YCE*YCE) 
CCCCCCCC        ALSHFT =  SQRT(XCE*XCE + YCE*YCE) 
          TETA=TETA1 + ATAN2(YCE,XCE)
          S=S1+ALSHFT*COS(TETA)
          Y=Y1+ALSHFT*SIN(TETA)
          CALL CALTRI(N,S,Y)
          TETA1 = TETA1 - PHI2RF

C         call fbgtxt
C         write(*,*) ' caltra mult ',xce,yce,ale,ATAN2(YCE,XCE),s1,s,y1,y

        ELSEIF(KPOS.EQ.2) THEN

        ELSE
c-------- KPOS=1

C         call fbgtxt
C         write(*,*) ' mult  S,Y,Z,TETA1,XI ',S,Y,Z,TETA1,XI

        ENDIF
C        DO II=4,10
C          IF(A(NOEL,II).NE.0.D0) AK=A(NOEL,II)/A(NOEL,3)/BORO *1.D4
C        ENDDO
C        DL=AL/2.D0

        PHI=0.D0
        CALL STORIG(NOEL,S,Y,Z,TETA1,XI,PHI)

        CALL SYNBOX(N,TETA1,S1,Y1,TETA,S,Y,AL,WIDTH)

        IF(KPOS.EQ.3) THEN
          XCE=0.D0
          YCE=0.D0
          ALE=A(NOEL,6+ND(NOEL)) 
          XIRF = 0.D0
          PHI2RF = -ALE 
          TETA=TETA1 
          S=S1
          Y=Y1
          CALL CALTRI(N,S,Y)
          TETA1 = TETA1 - PHI2RF

          XCE=0.D0   !!A(NOEL,4+ND(NOEL)) * CM2M
          YCE= A(NOEL,5+ND(NOEL)) * CM2M
          ALE= 0.D0
          XIRF = 0.D0
          PHI2RF = -ALE 
          ALSHFT = COEF * SQRT(XCE*XCE + YCE*YCE) 
          TETA=TETA1 + ATAN2(YCE,XCE)
          S=S1+ALSHFT*COS(TETA)
          Y=Y1+ALSHFT*SIN(TETA)
          CALL CALTRI(N,S,Y)
          TETA1 = TETA1 - PHI2RF

C         call fbgtxt
C         write(*,*) ' caltra mult ',xce,yce,ale,ATAN2(YCE,XCE),s1,s,y1,y

        ENDIF

      IF(N.EQ.1) THEN
        WRITE(6,FMT='(1X,A,1P,2G12.4,1X,2G12.4,1X,3G12.4)')
     >           TXT,SP,YP,S,Y,TETA,TETA1,WIDTH
      ELSEIF(N.EQ.2) THEN
C        CALL FBGTXT
C        WRITE(6,*) '  Now plotting : ',TXT
      ENDIF
      GOTO 998
C----- SEPARATEUR ELECTROSTATIQUE ANALYTIQUE
 35   CONTINUE
      IF(READAT) CALL RSEPAR
      IF(READAT) GOTO 998
C      CALL SEPARA(*99)
      GOTO 998
C----- RAZ LES COMPTEURS ( PEUT S'INTERCALER ENTRE PB CONSECUTIFS)
 36   CONTINUE
      GOTO 998
C----- IMPOSE  L'ORDRE  DE CALCUL  POUR LES ELMNTS KALC=3
C         ( ORDRE 2 PAR DEFAUT)
 37   CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1)
      IF(READAT) GOTO 998
      GOTO 998
C----- DIPOLE-M, DIPOLE ((i )2-D mid plane dipole map, (ii) analytical)
 39   CONTINUE
      IF(READAT) CALL RAIMAN(NDAT,NOEL,MXL,A,ND(NOEL))
      IF(READAT) GOTO 998
      GOTO 998
C----- RECHERCHE DU PLAN IMAGE ET DIMENSIONS D'IMAGE VERTICAUX
 41   CONTINUE
      GOTO 998
C----- RECHERCHE DU PLAN IMAGE ET DIMENSIONS D'IMAGE VERTICAUX
C----- FAISCEAU DISPERSE EN V
 42   CONTINUE
      GOTO 998
C----- DIMENSIONS ET POSITION VERTICALES DU FAISCEAU 
 43   CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1)
      IF(READAT) GOTO 998
      GOTO 998
C----- PLOTDATA PURPOSE STUFF
 44   CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1)
      IF(READAT) GOTO 998
      GOTO 998
C----- READS  THE  FORMATTED  FILE  'MYFILE' AND
C          COPIES  IT  INTO  THE  BINARY  FILE  'B_MYFILE',
C          OR RECIPROCAL.
 45   CONTINUE
      IF(READAT) CALL RBINAR
      IF(READAT) GOTO 998
      GOTO 998
C----- FIT
 102  CONTINUE
 46   CONTINUE
      CALL RFIT
      GOTO 998
C----- CARTE DE CHAMP CARTESIENNE MESUREE DU SPES3
C          D'APRES W. ROSCH, 1991
 47   CONTINUE
      KALC=2
      KUASEX=3
      IF(READAT) CALL RCARTE(1,2,ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- CARTE DE CHAMP CARTESIENNE 2-D DE CHALUT
 48   CONTINUE
      KALC=2
      KUASEX=4
      IF(READAT) CALL RCARTE(1,2,ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C-----  SPNTRK. Switch spin tracking
 49   CONTINUE
      IF(READAT) CALL RSPN
      IF(READAT) GOTO 998
      GOTO 998
C----- SPNPRT. PRINT SPIN STATES
 50   CONTINUE
      GOTO 998
C----- BEND. DIPOLE MAGNET
 51   CONTINUE
      IF(READAT) CALL RBEND(ND(NOEL))
      IF(READAT) GOTO 998
      XL = A(NOEL,10)
      AL = XL*COEF  * CM2M 
      REDUC=AL/4.D0
      WIDTH = D *REDUC

C      shift of origin (due to upstream fringe field extent) : 
      XI = -A(NOEL,20)  * CM2M  
      IF(A(NOEL,21).LT.0.D0) XI = 0.D0

      IF(A(NOEL,11).NE.0.D0) THEN
        WRITE(6,FMT='(/,A,/)')
     >    ' ** WARNING : cannot handle non zero skew angle in BEND'
        GOTO 998
      ENDIF
      B1  =A(NOEL,12)
      DR  =BORO/B1    * CM2M 
      KPOS = A(NOEL,70)
      TXT = 'BEND'
      A22 = A(NOEL,22)
      A42 = A(NOEL,42)
      A73 = A(NOEL,73)
      IF(KPOS.EQ.3) THEN
        IF(A73.NE.0.D0) THEN
          DEV2 = -A73
        ELSE
          DEV2 = ASIN(XL/2.D0/(BORO/B1))
        ENDIF
      ENDIF
      PHI2 = DEV2
      AL2=DR*SIN(PHI2)
     
      IF(KPOS.EQ.1) THEN
        COIN1 = +A22-TETREF
      ELSEIF(KPOS.EQ.2) THEN
        COIN1 = A73+A22-TETREF
      ELSEIF(KPOS.EQ.3) THEN
        COIN1 = A22-TETREF
      ENDIF

      TETREF = 0.D0
      NL=0
 510  CONTINUE
        NL=NL+1     
        IF(NOEL+NL.LE.NB) THEN
          IF(IQ(NOEL+NL).EQ.13) THEN          ! CHANGREF  
            TETREF=TETREF+A(NOEL+NL,3) * DEG2RD
            GOTO 510
          ENDIF
        ENDIF

      IF(KPOS.EQ.1) THEN
        COIN2 = +A42-TETREF
      ELSEIF(KPOS.EQ.2) THEN
        COIN2 = A73+A42-TETREF
      ELSEIF(KPOS.EQ.3) THEN
        COIN2 = A42-TETREF
      ENDIF
         
C       call fbgtxt
C         write(*,*) ' bend  S,Y,Z,TETA1,XI ',S,Y,Z,TETA1,XI

C      IDIP=1
C         CR1 = DR
C         CR2 = DR
         AL3=WIDTH*TAN(PHI2-COIN1)
         AL4=WIDTH*TAN(PHI2-COIN2)
         TETA1=TETA1-PHI2

         CALL STORIG(NOEL,S,Y,Z,TETA1-2.D0*PHI2,XI,2.D0*PHI2)

         CC=COS(TETA1)
         SS=SIN(TETA1)
         S0=S1+AL2*CC
         Y0=Y1+AL2*SS
         S=S0-(AL2+AL3)*CC-WIDTH*SS
         Y=Y0-(AL2+AL3)*SS+WIDTH*CC
         CALL CALTRI(N,S,Y)
         S=S0+(AL2+AL4)*CC-WIDTH*SS
         Y=Y0+(AL2+AL4)*SS+WIDTH*CC
         CALL CALTRI(N,S,Y)
         S=S0+(AL2-AL4)*CC+WIDTH*SS
         Y=Y0+(AL2-AL4)*SS-WIDTH*CC
         CALL CALTRI(N,S,Y)
         S=S0-(AL2-AL3)*CC+WIDTH*SS
         Y=Y0-(AL2-AL3)*SS-WIDTH*CC
         CALL CALTRI(N,S,Y)
         S=S0-AL2*CC
         Y=Y0-AL2*SS
         CALL CALTRI(N,S,Y)
         S=S0+AL2*CC
         Y=Y0+AL2*SS
         CALL CALTRI(N,S,Y)
         TETA1=TETA1-PHI2
         IF(N.EQ.1)  THEN
           WRITE(6,FMT='(1X,A,1P,2G12.4,1X,2G12.4,1X,3G12.4)')
     >            TXT,SP,YP,S,Y,TETA,TETA1,WIDTH

         ELSEIF(N.EQ.2) THEN
C           CALL TRTXT(S,Y,LABEL(NOEL,1),4,1)
C           CALL TRTXT(S0+AL2*CC/2.D0,Y0+AL2*SS/2.D0,TXT,4,1)

         ENDIF
      GOTO 998
C----- SOLENOIDE
 52   CONTINUE
      KALC=3
      KUASEX=20
      IF(READAT) CALL RSOLEN(ND(NOEL))
      IF(READAT) GOTO 998
      AL = A(NOEL,10) *COEF * CM2M 
      XI=0.D0
      PHI=0.D0
      CALL STORIG(NOEL,S1,Y1,Z,TETA1,XI,PHI)
      REDUC=AL/2.D0
      WIDTH = DS * REDUC  
      CALL SYNBOX(N,TETA1,S1,Y1,TETA,S,Y,AL,WIDTH)
      TXT = 'SOLE'
      IF(N.EQ.1) THEN
        WRITE(6,FMT='(1X,A,1P,2G12.4,1X,2G12.4,1X,3G12.4)')
     >           TXT,SP,YP,S,Y,TETA,TETA1,WIDTH
      ENDIF
C         call fbgtxt
C         write(*,*) ' SOLENO AL, S,Y,TETA1,width ',AL,S,Y,TETA1,width
      GOTO 998
C----- Plot transverse coordinates (for use on a workstation)
 53   CONTINUE
      GOTO 998
C----- PRINT STATE OF SPIN ON LOGICAL UNIT SPECIFIED
 54   CONTINUE
      IF(READAT) READ(NDAT,FMT='(A)') TA(NOEL,1)
      IF(READAT) GOTO 998
      GOTO 998
C----- CAVITE ACCELERATRICE
 55   CONTINUE
      IF(READAT) CALL RCAVIT
      IF(READAT) GOTO 998
C      CALL CAVITE(*99)
      GOTO 998
C----- DATA PARTICULE
 56   CONTINUE
C     .... MASSE(MeV/c2), CHARGE(C), G-yromagn., TPS VIE CDM(s), LIBRE
      IF(READAT) READ(NDAT,*) (A(NOEL,I),I=1,5)
      IF(READAT) GOTO 998
      GOTO 998
C----- COMMANDE DES SCALINGS ( B, FREQ...)
 57   CONTINUE
      IF(READAT) CALL RSCALE
      CALL FBGTXT
      WRITE(*,*) 
      WRITE(*,*) '-----------------------------------------------------'
      WRITE(*,*) '*** WARNING: Keyword "SCALING" appears in zgoubi.dat.'
      WRITE(*,*) ' Will not be accounted for. May entail mis-running,'
      WRITE(*,*) '   for instance wrong angles in BEND...'
      WRITE(*,*) '-----------------------------------------------------'
      WRITE(*,*) 
      I=IDLG('('' Press RETURN tp persue :'')','    ',1)
      IF(READAT) GOTO 998
      GOTO 998
C----- ELREVOL - CHAMP ELCTROSTATIQ Ex(R=0,X) 1-D MESURE. SYM CYLINDRIQ
 68   CONTINUE
      KFLD=LC
C----- BREVOL. CHAMP B(r=0,x) 1-D MESURE SUR AXE, A SYM CYLINDRIQ
 58   CONTINUE
      KALC=2
      KUASEX=8
      IF(READAT) CALL RCARTE(1,1,ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- POISSON. CARTE DE CHAMP CARTESIENNE 2-D FABRIQUEE PAR POISSON
 59   CONTINUE
      KALC=2
      KUASEX=5
      IF(READAT) CALL RCARTE(1,2,ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- DIPOLE. Similar to DIPOLES, but with a single magnet
 60   CONTINUE
      GOTO 39
C----- CARTE MESUREE SPECTRO KAON GSI (DANFISICS)
 61   CONTINUE
      KALC=2
      KUASEX=6
      IF(READAT) CALL RCARTE(1,2,ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- MAP2D: 2D B-FIELD MAP, NO SPECIAL SYMMETRY (P. Akishin, 07/1992)
 62   CONTINUE
      KALC=2
      KUASEX=9
      IF(READAT) CALL RCARTE(1,2,ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- B DECAPOLE ET DERIVEES CALCULES EN TOUT POINT (X,Y,Z)
 63   CONTINUE
      KALC=3
      KUASEX=5
      IF(READAT) CALL RQSOD(NDAT,NOEL,MXL,A,ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- B DODECAPOLE ET DERIVEES CALCULES EN TOUT POINT (X,Y,Z)
 64   CONTINUE
      KALC=3
      KUASEX=6
      IF(READAT) CALL RQSOD(NDAT,NOEL,MXL,A,ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- WIENFILT - INTEGR NUMERIQ.
 65   CONTINUE
      KALC=3
      KUASEX=21
      KFLD=ML
      IF(READAT) CALL RWIENF(ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- FAISTORE - Similar to FAISCNL, with additional options:
C        - Print every other run IPASS = mutltiple of IA
C        - Print after LABEL'ed elements
 66   CONTINUE
      IF(READAT) THEN
        READ(NDAT,*) TA(NOEL,1)
        READ(NDAT,*) A(NOEL,1)
      ENDIF
      GOTO 998
C----- SPNPRNLA - PRINT SPINS ON FILE, EACH IPASS=MULTIPL OF IA
 67   CONTINUE
      IF(READAT) THEN
        READ(NDAT,FMT='(A)') TA(NOEL,1)
        READ(NDAT,*) A(NOEL,1)
        IA=A(NOEL,1)
      ENDIF
      GOTO 998
C----- EL2TUB - LENTILLE ELECTROSTATIQ A 2 TUBES OU DIAPHRAG.
 69   CONTINUE
      KALC=3
      KUASEX=22
      KFLD=LC
      IF(READAT) CALL REL2TU(ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- UNIPOT. LENTILLE ELECTROSTATIQ A 3 TUBES OU DIAPHRAG.
 70   CONTINUE
      KALC=3
      KUASEX=23
      KFLD=LC
      IF(READAT) CALL RUNIPO(ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- ELMULT. E MULTIPOLAIRE ET DERIVEES CALCULES EN TOUT POINT (X,Y,Z)
 71   CONTINUE
      KALC=3
      KUASEX=MPOL+1
      KFLD=LC
      IF(READAT) CALL RMULTI(NDAT,NOEL,MXL,A,MPOL,ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- EBMULT. E+B MULTIPOLAIRES ET DERIVEES CALCULES EN TOUT POINT (X,Y,Z)
 72   CONTINUE
      KALC=3
      KUASEX=MPOL+1
      KFLD=ML
      IF(READAT) CALL REBMUL(ND(NOEL))
      IF(READAT) GOTO 998
C      CALL QUASEX(ND(NOEL))
      GOTO 998
 73   CONTINUE
C----- ESPACE LIBRE VIRTUEL
      IF(READAT) READ(NDAT,*) A(NOEL,1),A(NOEL,2)
      IF(READAT) GOTO 998
C      CALL ESLVIR(1,ZERO,ZERO)
      GOTO 998
C----- TRANSFORMATION FO(I,n)= A(,3)*FO(I,n) + A(,4)*FO(J,n) + A(,5)
 74   CONTINUE
      GOTO 998
C----- TRANSLATION  F(2,n)=F(2,n)+A(,1) et F(4,n)=F(4,n)+A(,2)
 75   CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1),A(NOEL,2)
      IF(READAT) GOTO 998
      GOTO 998
C----- POLARMES. CARTE DE CHAMP POLAIRE MESUREE
 76   CONTINUE
      KALC=2
      KUASEX=22
      IF(READAT) CALL RCARTE(2,2,ND(NOEL))
      IF(READAT) GOTO 998
C      CALL AIMANT(ND(NOEL))
      GOTO 998
C----- TRANSLATION X,Y,Z ET ROTATION RX,RY,RZ
 77   CONTINUE
      IF(READAT) READ(NDAT,*) (A(NOEL,II),II=1,6)
      IF(READAT) GOTO 998
C      CALL TRAROT(A(NOEL,1),A(NOEL,2),A(NOEL,3)
C     >  ,A(NOEL,4)*RAD,A(NOEL,5)*RAD,A(NOEL,6)*RAD)
      GOTO 998
C----- SRLOSS. Switch synchrotron ligth on
 78   CONTINUE
      IF(READAT) CALL RSRLOS
      GOTO 998
C----- AVERAGEORB. Calculate average orbit from multiturn 
 79   CONTINUE
      IF(READAT) CALL RAVORB
      GOTO 998
C----- OPTICS. Calculate Twiss functns, transport and print at
 80   CONTINUE
      GOTO 998
C-----  GASCAT. Switch gas-scattering
 81   CONTINUE
      GOTO 998
C----- UNDULATOR
 82   CONTINUE
      KALC=1
      KUASEX=9
      IF(READAT) CALL RUNDUL(ND(NOEL))
      GOTO 998
C----- ELCYLDEF
 83   CONTINUE
      STOP ' Element ELCYLDEF not implemented'
C      IF(READAT) CALL RELCYL(ND(NOEL))
C      IF(READAT) GOTO 998
C      KALC = 3
C      KUASEX = 24
C      KFLD=LC
      GOTO 998
C----- ELMIR
 84   CONTINUE
      IF(READAT) CALL RELMIR(ND(NOEL))
      IF(READAT) GOTO 998
      IF(A(NOEL,10+INT(A(NOEL,10))+2).LE.11) THEN
C------- Mirror
        XL = A(NOEL,42)
        AL = XL*COEF  * CM2M 
        REDUC=AL/2.D0
        WIDTH = D *REDUC
      ELSE
C------- Transmission lens
        XL = 0.D0
        NPLATE=INT(A(NOEL,10))
        DO 840 I=1,NPLATE
 840      XL = XL + A(NOEL,10+I)
        AL = XL*COEF  
        REDUC=AL
        IF(A(NOEL,10+INT(A(NOEL,10))+2).EQ.22) THEN
C------- LH
          WIDTH = .5D0*A(NOEL,11+NPLATE)
        ELSE
C------- LV
          WIDTH = DQ*REDUC
        ENDIF
C        AK = A(NOEL,21) - A(NOEL,20)
        TXT = 'LENS'
        GOTO 310
      ENDIF
      XI = -A(NOEL,11)
      IF(A(NOEL,43).EQ.0.D0) THEN
        WRITE(6,FMT='(/,A,/)')
     >    ' ** WARNING : cannot handle zero ALE in ELMIR'
        GOTO 998
      ENDIF
      TXT = 'ELMI'
      KPOS=A(NOEL,40)
      IF(KPOS.EQ.3) THEN
CCCCCCCCCCC        DEV2 = -A(NOEL,43)
        DEV2 = A(NOEL,43)
      ELSE
        WRITE(6,FMT='(/,A,/)')
     >    ' ** WARNING : USE KPOS=3 in ELMIR mirror'
        GOTO 998
      ENDIF

      TETREF = 0.D0
        COIN1 = DEV2
        COIN2 = COIN1
         
      PHI2 = DEV2

        IF(DEV2.LT.0.D0) WIDTH=-WIDTH
         AL2=AL/2.D0
         AL3=WIDTH*TAN(PHI2-COIN1)
         AL4=WIDTH*TAN(PHI2-COIN2)
         TETA1=TETA1-PHI2

      CALL STORIG(NOEL,S,Y,Z,TETA1,XI,2.D0*PHI2)

         CC=COS(TETA1)
         SS=SIN(TETA1)
         S0=S1+AL2*CC
         Y0=Y1+AL2*SS

         S=S0-(AL2+AL3)*CC-WIDTH*SS*2.D0
         Y=Y0-(AL2+AL3)*SS+WIDTH*CC*2.D0
         CALL CALTRI(N,S,Y)
         S=S0+(AL2+AL4)*CC-WIDTH*SS*2.D0
         Y=Y0+(AL2+AL4)*SS+WIDTH*CC*2.D0
         CALL CALTRI(N,S,Y)
         S=S0+AL2*CC
         Y=Y0+AL2*SS
         CALL CALTRI(N,S,Y)
         S=S0-AL2*CC
         Y=Y0-AL2*SS
         CALL CALTRI(N,S,Y)
         S=S0-(AL2+AL3)*CC-WIDTH*SS
         Y=Y0-(AL2+AL3)*SS+WIDTH*CC
         CALL CALTRI(N,S,Y)
         S=S0+(AL2+AL4)*CC-WIDTH*SS
         Y=Y0+(AL2+AL4)*SS+WIDTH*CC
         CALL CALTRI(N,S,Y)
         S=S0+AL2*CC
         Y=Y0+AL2*SS
         CALL CALTRI(N,S,Y)

         TETA1=TETA1-PHI2
         IF(N.EQ.1)  THEN
         WRITE(6,FMT='(1X,A,1P,2G12.4,1X,2G12.4,1X,3G12.4)')
     >            TXT,SP,YP,S,Y,TETA,TETA1,WIDTH
         ELSEIF(N.EQ.2) THEN
C           CALL TRTXT(S,Y,LABEL(NOEL,1),4,1)
C           CALL TRTXT(S0+AL2*CC/2.D0,Y0+AL2*SS/2.D0,TXT,4,1)
         ENDIF
      GOTO 998
C----- ELCMIR
 85   CONTINUE
      KALC=3
      KUASEX=26
      KFLD=LC
      IF(READAT) CALL RELCMI(ND(NOEL))
      GOTO 998
C----- MAP2D_E: 2D E-FIELD MAP
 86   CONTINUE
      KFLD=LC
      GOTO 62
C----- SRPRNT. Print/Store S.R. loss tracking statistics into file
 87   CONTINUE
      GOTO 998
C----- BETATRON. Betatron core
 88   CONTINUE
      IF(READAT) READ(NDAT,*) DPKCK
      GOTO 998
C----- TWISS. Compute linear lattice functions, chromaticity, etc. 
 89   CONTINUE
C                                        Fac_dp   Fac-ampl
      IF(READAT) READ(NDAT,*) A(NOEL,1),A(NOEL,2),A(NOEL,3)
      GOTO 998
C----- END. End of run, except for some options that may need more
 90   CONTINUE
      GOTO 995
C----- FFAG. FFAG Sector multi-dipole. 
 91   CONTINUE
      KALC = 1
      KUASEX = 27
      IF(READAT) THEN
        CALL RFFAG(ND(NOEL))
      ELSE
      ENDIF
      GOTO 998
C----- HELIX. Helical field (twisted dipole)
 92   CONTINUE
      KALC = 3
      KUASEX = 28
      IF(READAT) THEN
        CALL RHELIX(ND(NOEL))
      ELSE
      ENDIF
      GOTO 998
C-----  CSR. Switch  coherent SR interaction
 93   CONTINUE
      IF(READAT) CALL RCSR
      GOTO 998
C----- Each particle is pushed a Delta-S distance
 94   CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1)
      GOTO 998
C----- COILS. 
 95   CONTINUE
      KALC = 3
      KUASEX=29
      IF(READAT) THEN 
        CALL RCOILS(NDAT,NOEL,MXL,A,ND(NOEL))
      ELSE
      ENDIF
      GOTO 998
C----- GETFITVAL.  Get parameter values resulting from FIT, stored in zgoubi.fitVal
 96   CONTINUE
      IF(READAT) THEN
        CALL RFITGT
      ENDIF
      GOTO 998
C----- SUPERPOSE.  To superimpose magnets. 2B developped
 97   CONTINUE
C      IF(READAT) CALL RSUPER(NDAT,NOEL,MXL,A)
      KUASEX = 31
      GOTO 998
C----- MARKER. 
 98   CONTINUE
      IF(LABEL(NOEL,2).EQ.'.plt') THEN
        AL = 0.D0
        XI=0.D0
        PHI=0.D0
        CALL STORIG(NOEL,S,Y,Z,TETA1,XI,PHI)
        TXT = 'MARKER  / LABEL2=.plt'
        IF(N.EQ.1) 
     >    WRITE(6,FMT='(1X,A,1P,2G12.4,1X,2G12.4,1X,3G12.4)')
     >           TXT,SP,YP,S,Y,TETA,TETA1,WIDTH
      ENDIF
      GOTO 998
C----- DIPOLES. A set of neiboring or overlapping dipoles. 
 99   CONTINUE
      KALC =1
      KUASEX=32
      IF(READAT) CALL RDIPS(ND(NOEL))
      GOTO 998
C----- TRACKING. 
 100  CONTINUE
        READ(NDAT,*) NLA, NLB 
      GOTO 998
C----- FFAG-SPI. FFAG, spiral. 
 101  CONTINUE
      IF(READAT) CALL RFFAGS(ND(NOEL),NDAT,
     >    NMAG,AT,RM,ACN,OP,QSIE,OM,QSIS,RE,TE,RS,TS)
       mod = 0
      CALL READCC(MOD,RM*CM2M,RM*CM2M)
      IF(READAT) GOTO 998

      if(n.eq.2) 
     >CALL PRFFAG(NMAG,AT,RM,ACN,OP,QSIE,OM,QSIS,RE,TE,RS,TS)

      XL = 2.d0 * RM * sin((op - om)*deg2rd/2.d0)  
      AL = XL*COEF  * CM2M 
      REDUC=AL/4.D0
      WIDTH = D *REDUC

C      shift of origin (due to upstream fringe field extent) : 
      AMIN = AT/2.D0 +ACN       
      XI =  RM * TAN(op-(2.d0*acn-at/2.d0)*DEG2RD) * CM2M
      TXT = 'FFAG-SPI'

      B1  =A(NOEL,12)
      DR  = rm    * CM2M 
          DEV2 = 0.5d0 * (at*deg2rd + te - ts)
C           write(*,*) at,te,ts,dev2/deg2rd
      PHI2 = DEV2
      AL2=DR*SIN(PHI2)
     
      Adrift = tetF - AF
      COIN1 = Adrift/2.D0 - QSie
      COIN2 = Adrift/2.D0 + QSis

      TETREF = 0.D0

         AL3=WIDTH*TAN(PHI2-COIN1)
         AL4=WIDTH*TAN(PHI2-COIN2)
         TETA1=TETA1-PHI2

         CALL STORIG(NOEL,S,Y,Z,TETA1-2.D0*PHI2,XI,2.D0*PHI2)

         CC=COS(TETA1)
         SS=SIN(TETA1)
         S0=S1+AL2*CC
         Y0=Y1+AL2*SS
         S=S0-(AL2+AL3)*CC-WIDTH*SS
         Y=Y0-(AL2+AL3)*SS+WIDTH*CC
         CALL CALTRI(N,S,Y)
         S=S0+(AL2+AL4)*CC-WIDTH*SS
         Y=Y0+(AL2+AL4)*SS+WIDTH*CC
         CALL CALTRI(N,S,Y)
         S=S0+(AL2-AL4)*CC+WIDTH*SS
         Y=Y0+(AL2-AL4)*SS-WIDTH*CC
         CALL CALTRI(N,S,Y)
         S=S0-(AL2-AL3)*CC+WIDTH*SS
         Y=Y0-(AL2-AL3)*SS-WIDTH*CC
         CALL CALTRI(N,S,Y)
         S=S0-AL2*CC
         Y=Y0-AL2*SS
         CALL CALTRI(N,S,Y)
         S=S0+AL2*CC
         Y=Y0+AL2*SS
         CALL CALTRI(N,S,Y)
         TETA1=TETA1-PHI2
         IF(N.EQ.1)  THEN
           WRITE(6,FMT='(1X,A,1P,2G12.4,1X,2G12.4,1X,3G12.4)')
     >            TXT,SP,YP,S,Y,TETA,TETA1,WIDTH

         ELSEIF(N.EQ.2) THEN
C           CALL TRTXT(S,Y,LABEL(NOEL,1),4,1)
C           CALL TRTXT(S0+AL2*CC/2.D0,Y0+AL2*SS/2.D0,TXT,4,1)

         ENDIF
      GOTO 998
C----- EMMA. Read  2-D or 3-D field map (e.g., as obtained from TOSCA code), 
C      with mesh either cartesian (KART=1) or cylindrical (KART=2). 
 103  CONTINUE
      GOTO 998

c      MOD = 33
c      CALL READCC(MOD,RM*CM2M,RM*CM2M)
c      CALL PRFFAG(NOEL,NMAG,AT,RM,ACN,OP,QSIE,OM,QSIS,RE,TE,RS,TS)

      GOTO 998
      END
