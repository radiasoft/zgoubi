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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE ZGOUBI(NL1,NL2,READAT,
     >                                 NBLMN,ENDFIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL READAT, ENDFIT

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     >     IREP(MXT),AMQLU,PABSLU
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LABEL
      COMMON /LABEL/ LABEL(MXL,2)
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE 'MXFS.H'
      INCLUDE 'MXSCL.H'     
      COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLEY
      SAVE KLEY
      CHARACTER(80) TITRE
      COMMON/TITR/ TITRE 

      PARAMETER(MPOL=10)
      DIMENSION ND(MXL)

C----- For printing after occurence of pre-defined labels
      PARAMETER(MLB=10)
      CHARACTER(LBLSIZ) LBL(MLB), LBLSP(MLB)
      LOGICAL PRLB, PRLBSP
      SAVE KPRT, PRLB, KPRTSP, PRLBSP

C----- Pick-up signal
      PARAMETER (MXPUD=9,MXPU=5000)
      COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      PARAMETER (MPULAB=5)
      CHARACTER(10) PULAB
      COMMON/COT/ PULAB(MPULAB)

C----- Set to true by REBELOTE if last turn to be stopped at NOELB<MAX_NOEL
      LOGICAL REBFLG
      SAVE REBFLG, NOELRB

C----- Tells whether FIT is active or not
      LOGICAL FITING
C----- To get values into A(), from earlier FIT
      LOGICAL FITGET
      SAVE FITGET

      LOGICAL TOMANY, STRACO
      CHARACTER(80) TXTEMP
      CHARACTER(1) TXT1
      INTEGER DEBSTR,FINSTR

      CHARACTER(80) TXTELT, TXTELO
      SAVE TXTELT

      CHARACTER(300) SYSCMD

      INCLUDE 'PARIZ.H'
      INCLUDE 'FILPLT.H'

      PARAMETER (I0=0, I1=1, I2=2, I3=3, I5=5, I6=6)

      CHARACTER(LBLSIZ) LBLOPT
      SAVE KOPIMP, KOPTCS, LBLOPT
      LOGICAL OKLNO

      LOGICAL EMPTY, IDLUNI, OKCPLD

      SAVE LNOPTI
      SAVE PNLTGT
      SAVE OKLNO

C This INCLUDE must stay located right before the first statement
      CHARACTER(KSIZ) KLEO
 
      LOGICAL PRDIC

      LOGICAL FITFNL
      SAVE FITFNL

      CHARACTER(132) TXT132
      LOGICAL STRCON 

      DATA PRDIC / .FALSE. /
      DATA OKLNO / .FALSE. /

      INCLUDE 'LSTKEY.H'
      
      IF(ENDFIT) GOTO 998 

      IF(NL2 .GT. MXL) CALL ENDJOB(
     >      'Too  many  elements  in  the  structure, max is',MXL)

      IF(READAT) THEN
        CALL PRDATA(
     >              LABEL,NBLMN)
        CALL FITSTA(I5,
     >                 FITING)
        CALL FITST2(NBLMN)
        IF(.NOT.FITING) NL2=NBLMN
        CALL LINGUA(LNG)
        CALL RESET
C------- Print after defined labels. Switched on by FAISTORE.
        PRLB = .FALSE.
C------- Print after defined labels. Switched on by SPNSTORE.
        PRLBSP = .FALSE.
      ENDIF
        
C----- Get FIT status
      CALL FITSTA(I5,
     >               FITING)
      IF(FITING) CALL RESET2

      IF(READAT) READ(NDAT,503) TITRE
 503  FORMAT(A)

CCCCCCCCCCCCfor LHC : do    REWIND(4)

      IF(.NOT. FITING) THEN 
        IF(NRES .GT. 0) WRITE(6,905) TITRE
 905  FORMAT(/,1X,'Title : ',A80,//)
      ENDIF

      TOMANY = .FALSE.
      NOEL = NL1-1

 998  CONTINUE
  
C YD FM. 28/01/2014
      IF(NOEL .GT. 0) THEN

       IF(PRLB) THEN
C------- Print after Lmnt with defined LABEL - from Keyword FAISTORE
C        LBL contains the LABEL['s] after which print shall occur
        IF( STRACO(NLB,LBL,LABEL(NOEL,1),
     >                                   IL) 
     >    .OR. LBL(1).EQ.'all' .OR. LBL(1).EQ.'ALL') 
     >    CALL IMPFAI(KPRT,NOEL,KLE(IQ(NOEL)),LABEL(NOEL,1),
     >                                              LABEL(NOEL,2)) 
       ENDIF
       IF(PRLBSP) THEN
C------- Print after Lmnt with defined LABEL - from Keyword SPNSTORE
C        LBLSP contains the LABEL['s] after which print shall occur
        IF( STRACO(NLB,LBLSP,LABEL(NOEL,1),
     >                                     IL) ) 
     >    CALL SPNPRN(KPRTSP,NOEL,KLE(IQ(NOEL)),LABEL(NOEL,1),
     >                                              LABEL(NOEL,2)) 
       ENDIF

       IF(KCO .EQ. 1) THEN
C------- Calculate pick-up signal
C        PULAB contains the NPU LABEL's at which CO is calculated 
        IF( STRACO(NPU,PULAB,LABEL(NOEL,1),
     >                                  IL) ) 
     >    CALL PCKUP(NOEL)
C     >    CALL PCKUP(NOEL,KLE(IQ(NOEL)),LABEL(NOEL,1),LABEL(NOEL,2))

       ENDIF


       IF(KOPTCS .EQ. 1) THEN
C------- Transport beam matrix and print at element ends. Switched by OPTICS keyword.

        IF(KUASEX .NE. -99) THEN
C         Only if keyword is of optical element type or [MC]OBJET or END or MARKER

          IF(EMPTY(LBLOPT) .OR. 
     >      LBLOPT .EQ. 'ALL' .OR. LBLOPT .EQ. 'all' .OR. 
     >              LBLOPT.EQ.LABEL(NOEL,1)) THEN
C            CALL OPTICC(LNOPTI,NOEL,KOPIMP,PRDIC,OKCPLD)
            CALL OPTICC(        NOEL,KOPIMP,PRDIC,OKCPLD)

          ENDIF
        ENDIF
       ENDIF

C This was introduced so to restrain OPTICS to optics-type keywords.  Any additional desired 
C location should be assigned KUASEX=0 for OPTICC to operate (e.g., DRIFT, MARKER...)
          KUASEX = -99


       IF(REBFLG) THEN
C----- Set to true by REBELOTE : last turn to be stopped at NOELB<MAX_NOEL
        IF(IPASS.EQ.NRBLT+1) THEN
          CALL REBEL7(
     >                NOELB)
          IF(NOEL.EQ.NOELB) THEN
            IKLE = 11  
            KLEY = KLE(IKLE)   ! 'REBELOTE'
            NOEL = NOELRB
            IQ(NOEL) = IKLE
            GOTO 187
          ENDIF
        ENDIF
       ENDIF

      ENDIF ! NOEL.GT.0 

      IF(READAT) THEN
 188    READ(NDAT,*,ERR=999,END=999) KLEY
        IF(KLEY(DEBSTR(KLEY):DEBSTR(KLEY)) .EQ. '!') GOTO 188
        DO IKLE=1,MXKLE
          IF(KLEY .EQ. KLE(IKLE)) THEN
            NOEL = NOEL+1
            IF( NOEL .EQ. MXL+1) THEN
              TOMANY=.TRUE.
              IF(NRES .GT. 0) THEN
                WRITE(NRES,*) ' PROCEDURE STOPPED: too many elements'
                WRITE(NRES,*) ' (number of elements should not exceed '
     >                                                       ,MXL,').'
                WRITE(NRES,*) ' Increase  MXL  in   MXLD.H'
              ENDIF
              CALL ENDJOB(' Increase  MXL  in   MXLD.H',-99)
            ENDIF
            IQ(NOEL) =  IKLE
            GOTO 187
          ENDIF
        ENDDO
        GOTO 999
      ELSE
C------- Gets in case of "FIT"
        IF (NOEL .EQ. NL2 ) RETURN
        NOEL = NOEL+1
        IKLE = IQ(NOEL)
        KLEY = KLE(IKLE)
      ENDIF
 
 187  CONTINUE
      
      IF(NRES.GT.0) THEN
        WRITE(NRES,201)
 201    FORMAT(/,132('*'))
        WRITE(NRES,334) NOEL,'Keyword, label(s) :',
     >  KLEY,LABEL(NOEL,1),LABEL(NOEL,2)
 334    FORMAT(2X,I5,4(2X,A),/)
        CALL FLUSH2(NRES,.FALSE.)
        WRITE(TXTELT,FMT='(I5,A1,I5,1X,A10,2(A1,A))') 
     >    NOEL,'/',NBLMN,KLEY,'/',LABEL(NOEL,1),'/',LABEL(NOEL,2)
        IF(IPASS.EQ.1) CALL ARRIER(TXTELT)
      ENDIF
 
      KFLD=MG

CCCCCCCCCCCCfor LHC ; do REWIND(4).
C      IF(1000*(NOEL/1000) .EQ. NOEL) REWIND(4)
C            write(88,*) ' zgoubi ',NDAT,NRES,NPLT,NFAI,NMAP,NSPN,ipass

C Go to "GOTO(...) IKLE" ---------------------------
      GOTO 1001
  999 CONTINUE
C---------------------------------------------------

      IF(NRES.GT.0) THEN
        WRITE(NRES,201)

        WRITE(NRES,200) KLEY
 200    FORMAT(/,10X,' MAIN PROGRAM : Execution ended upon key  ',A)
        WRITE(NRES,201) 
      ENDIF

      RETURN 
 
C----- DRIFT, ESL. Free space. 
 1    CONTINUE
      KUASEX = 0
      IF(READAT) READ(NDAT,*) A(NOEL,1)
      IF(FITGET) CALL FITGT1
      CALL ESL(1,A(NOEL,1),1,IMAX)
      GOTO 998
C----- AIMANT. Dipole with computed field map in cylindrical coordinates
 2    CONTINUE
      IF(READAT) CALL RAIMAN(NDAT,NOEL,MXL,A,ND(NOEL))
      IF(FITGET) CALL FITGT1
      KALC = 2
      KUASEX = 20
      CALL AIMANT(ND(NOEL))
      GOTO 998
C----- QUADRUPO - B quadrupolaire et derivees calcules en tout point (X,Y,Z)
 3    CONTINUE
      KALC = 3
      KUASEX = 2
      IF(READAT) THEN
        CALL RQSOD(NDAT,NOEL,MXL,A,ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- SEXTUPOL - B SEXTUPOLAIRE  ET DERIVEES CALCULES EN TOUT POINT (X,Y,Z)
4     CONTINUE
      KALC = 3
      KUASEX = 3
      IF(READAT) THEN
        CALL RQSOD(NDAT,NOEL,MXL,A,ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- IMAGE. RECHERCHE DU PLAN IMAGE ET DIMENSIONS D'IMAGE HORIZONTAUX
C----- FAISCEAU MONOCHROMATIQUE
 5    CALL FOCALE(1)
      GOTO 998
C----- IMAGES. RECHERCHE DU PLAN IMAGE ET DIMENSIONS D'IMAGE HORIZONTAUX
C----- FAISCEAU POLYCHROME
 6    CALL FOCALE(2)
      GOTO 998
C----- FAISCEAU. Print current beam in zgoubi.res
7     IF(NRES.GT.0) CALL IMPTRA(1,IMAX,NRES)
      GOTO 998
C----- FAISCNL. Stores beam at current position (in .fai type file)
 8    CONTINUE
      IF(READAT) CALL RFAIST(I0,
     >                         PRLB,KPRT,LBL,NLB)
      IF(TA(NOEL,1).NE.'none') THEN
        NLB = 0
        TXTEMP = TA(NOEL,1)
        TXTEMP=TXTEMP(DEBSTR(TXTEMP):FINSTR(TXTEMP))
        CALL IMPFAW(TXTEMP,LBL,NLB)
        CALL IMPFAI(I0,NOEL-1,KLE(IQ(NOEL-1)),LABEL(NOEL-1,1),
     >                                             LABEL(NOEL-1,2))
      ENDIF
      GOTO 998
C----- TRAVERSEE D'UNE CIBLE
9     CONTINUE
      IF(READAT) CALL RCIBLE
      IF(FITGET) CALL FITGT1
      CALL CIBLE
      GOTO 998
C----- FOCALE. DIMENSIONS DU FAISCEAU @ XI
10    CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1)
      IF(FITGET) CALL FITGT1
      CALL FOCALE(3)
      GOTO 998
C----- REBELOTE. Passes NRBLT more times thru the structure
11    CONTINUE
      IF(READAT) CALL RREBEL(LABEL,kle)
      CALL REBEL(READAT,KLE,LABEL,
     >                            REBFLG,NOELRB)
      CALL KSMAP0
      ENDFIT = .FALSE.   ! Used in zgoubi_main. Purpose : make REBELOTE compatible with FIT
      GOTO 998
C----- QUADISEX. Champ creneau B = B0(1+N.Y+B.Y2+G.Y3) plan median
 12   CONTINUE
      KALC =1
      KUASEX = 3
      IF(READAT) CALL RSIMB(ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL  QUASEX(ND(NOEL))
      GOTO 998
C----- CHANGREF. Translation X,Y et rotation du referentiel courant
 13   CONTINUE
      IF(READAT) CALL RCHANG
      IF(FITGET) CALL FITGT1
      CALL CHREFE
      GOTO 998
C----- SEXQUAD. Champ CRENEAU B = B0(N.Y+B.Y2+G.Y3) PLAN MEDIAN
 14   CONTINUE
      KALC =1
      KUASEX = 4
      IF(READAT) CALL RSIMB(ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL  QUASEX(ND(NOEL))
      GOTO 998
C----- AIMANT TOROIDAL
 15   CONTINUE
      KALC =1
      KUASEX = 7
      IF(READAT) CALL RSIMB(ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL  QUASEX(ND(NOEL))
      GOTO 998
C----- OBJETA. MONTE-CARLO GENERATION OF MU OR PI FROM ETA DECAY 
C          ( D'APRES Benjamin MAYER , FEVRIER 1990 )
 17   CONTINUE
      IF(READAT) CALL ROBJTA
      CALL OBJETA
      GOTO 998
C----- MATRIX. COEFFICIENTS D'ABERRATION A L'ABSCISSE COURANTE
 18   CONTINUE
      IF(READAT) CALL RMATRX
      IF(FITGET) CALL FITGT1
      OKCPLD = NINT(A(NOEL,4)) .EQ. 1
      CALL MATRIC(NINT(
     >  A(NOEL,1)),NINT(A(NOEL,2)),NINT(A(NOEL,3)),OKCPLD)
      GOTO 998
C----- CHAMBR. Stops and records trajectories out of chamber limits
 19   CONTINUE
      IF(READAT) CALL RCOLLI
      IF(FITGET) CALL FITGT1
      CALL CHAMB
      GOTO 998
C----- Champ Q-POLE SPECIAL PROJET SPES2
 20   CONTINUE
      KALC =1
      KUASEX = 1
C ............... ADD   IF(READAT) CALL ............
      IF(FITGET) CALL FITGT1
      CALL  QUASEX(ND(NOEL))
      GOTO 998
C----- CARTEMES. CARTE DE Champ CARTESIENNE MESUREE DU SPES2
 21   CONTINUE
      KALC =2
      KUASEX = 1
      IF(READAT) CALL RCARTE(I1,I2,
     >                             ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- CHANGE LE FAISCEAU Y,T,Z,P EN -Y,-T,-Z,-P
C       ( EQUIVALENT A CHANGEMENT DE SIGNE DU Champ DIPOLAIRE }
 22   CONTINUE
      CALL YMOINY
      GOTO 998
C----- B OCTUPOLAIRE ET DERIVEES CALCULES EN TOUT POINT (X,Y,Z)
 23   CONTINUE
      KALC = 3
      KUASEX = 4
      IF(READAT) THEN
        CALL RQSOD(NDAT,NOEL,MXL,A,ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- OBJET. 
24    CONTINUE
      KUASEX = 0
      IF(READAT) CALL ROBJET
      IF(FITGET) CALL FITGT1
      CALL OBJETS
      GOTO 998
C----- MCOBJET. Object defined by Monte-Carlo
 25   CONTINUE
      KUASEX = 0
      IF(READAT) CALL RMCOBJ
      IF(FITGET) CALL FITGT1
C FM July 2014 - Inhibited call to mcobja for eRHIC, to allow varying Z0,P0. 
C      IF(.NOT. FITING) THEN
        CALL MCOBJ
C      ELSE
C        CALL MCOBJA
C        CALL CNTRST
C      ENDIF
      GOTO 998
C----- DESINTEGRATION EN COURS DE VOL
26    CONTINUE
      IF(READAT) CALL RMCDES
      IF(.NOT. FITING)  CALL MCDESI(*999)
      IF(FITGET) CALL FITGT1
      GOTO 998
C----- HISTOGRAMME DE Y T Z P D
27    CONTINUE
      IF(READAT) CALL RHIST(NDAT,NOEL,MXL,A,TA)
      IF(FITGET) CALL FITGT1
      CALL HIST
      GOTO 998
C----- TRANSMAT. TRANSFERT MATRICIEL AU SECOND ORDRE
 28   CONTINUE
      IF(READAT) CALL RTRANS
      IF(FITGET) CALL FITGT1
      CALL TRANSM
      GOTO 998
C----- AIMANT VENUS (Champ CONSTANT DANS UN RECTANBLE)
 29   CONTINUE
      KALC =1
      KUASEX = 5
      IF(READAT) CALL RSIMB(ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL  QUASEX(ND(NOEL))
      GOTO 998
C----- AIMANT PS170 (Champ CONSTANT DANS UN CERCLE)
 30   CONTINUE
      KALC =1
      KUASEX = 6
      IF(READAT) CALL RSIMB(ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL  QUASEX(ND(NOEL))
      GOTO 998
C----- TOSCA. Read  2-D or 3-D field map (e.g., as obtained from TOSCA code), 
C      with mesh either cartesian (KART=1) or cylindrical (KART=2). 
 31   CONTINUE
      KALC = 2
      IF(READAT) CALL RCARTE(KART,I3,
     >                               ND(NOEL))
      IF    (A(NOEL,22) .EQ. 1) THEN
C        KZMA = 1 ; 2-D map
        KUASEX = 2
      ELSEIF(A(NOEL,22) .GT. 1) THEN
C        KZMA > 1 ; 3-D map
        KUASEX = 7
        IF(IZ.LE.1) CALL ENDJOB(' *** ERROR ; cannot use a 3-D map, need
     >  recompile zgoubi, using IZ>1 in PARIZ.H',-99) 
      ENDIF
      IF(FITGET) CALL FITGT1
      IF    (KART.EQ.1) THEN 
        CALL QUASEX(ND(NOEL))
      ELSEIF(KART.EQ.2) THEN 
        CALL AIMANT(ND(NOEL))
      ENDIF
      GOTO 998
C----- AUTOREF. CHGNT RFRNTIEL AUTMTIQ => PLAN NORMAL A LA TRAJ. DE REFERENCE
 32   CONTINUE
      IF(READAT) CALL RAUTOR
      IF(FITGET) CALL FITGT1
      CALL AUTORF
      GOTO 998
C----- COLLIMA. Collimator ( particules out => IEX = -4 )
 33   CONTINUE
      IF(READAT) CALL RCOLLI
      IF(FITGET) CALL FITGT1
      CALL COLLIM
      GOTO 998
C----- MULTIPOL. B ET DERIVEES CALCULES EN TOUT POINT (X,Y,Z)
 34   CONTINUE
      KALC = 3
      KUASEX = MPOL+1
      IF(READAT) THEN
        CALL RMULTI(NDAT,NOEL,MXL,A,MPOL,ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- SEPARATEUR ELECTROSTATIQUE ANALYTIQUE
 35   CONTINUE
      IF(READAT) CALL RSEPAR
      IF(FITGET) CALL FITGT1
      CALL SEPARA(*999)
      GOTO 998
C----- RAZ LES COMPTEURS ( PEUT S'INTERCALER ENTRE PB CONSECUTIFS)
 36   CONTINUE
      CALL RESET
      GOTO 998
C----- IMPOSE  L'ORDRE  DE CALCUL  POUR LES ELMNTS KALC = 3
C         ( ORDRE 2 PAR DEFAUT)
 37   CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1)
      IF(FITGET) CALL FITGT1
      CALL MODORD
      GOTO 998
C----- DIPOLE-M. 2-D mid plane fabricated dipole map. Polar coordinates.
 39   CONTINUE
      IF(READAT) CALL RAIMAN(NDAT,NOEL,MXL,A,ND(NOEL))
      KALC = 2
      KUASEX = 21
      IF(FITGET) CALL FITGT1
      CALL AIMANT(ND(NOEL))
      GOTO 998
C----- IMAGEZ. RECHERCHE DU PLAN IMAGE ET DIMENSIONS D'IMAGE VERTICAUX
 41   CALL FOCALE(-1)
      GOTO 998
C----- IMAGESZ. RECHERCHE DU PLAN IMAGE ET DIMENSIONS D'IMAGE VERTICAUX
C----- FAISCEAU DISPERSE EN V
 42   CALL FOCALE(-2)
      GOTO 998
C----- FOCALEZ. DIMENSIONS ET POSITION VERTICALES DU FAISCEAU @ XI
 43   CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1)
      IF(FITGET) CALL FITGT1
      CALL FOCALE(-3)
      GOTO 998
C----- PLOTDATA PURPOSE STUFF
 44   CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1)
      IF(FITGET) CALL FITGT1
      CALL PLTDAT
      GOTO 998
C----- READS  THE  FORMATTED  FILE  'MYFILE' AND
C          COPIES  IT  INTO  THE  BINARY  FILE  'B_MYFILE',
C          OR RECIPROCAL.
 45   CONTINUE
      IF(READAT) CALL RBINAR
      CALL BINARY
      GOTO 998
C----- FIT. FIT2. Two methods are available 
 102  MTHOD = 2
      GOTO 461
 46   CONTINUE
      MTHOD = 1
 461  CONTINUE
      CALL FITNU2(MTHOD)
      IF(READAT) CALL RFIT(KLEY,
     >                         PNLTGT,ICPTMA,FITFNL)
      CALL FITNU4(FITFNL)
      IF(MTHOD.EQ.1) CALL MINO12(PNLTGT,ICPTMA)
      IF(MTHOD.EQ.2) CALL NMMIN2(PNLTGT,ICPTMA)
      CALL CPTFC1(ICPTMA)
      FITING = .TRUE.
      CALL FITSTA(I6,FITING)
      CALL FITST2(NOEL)
      FITGET = .FALSE.
      RETURN
C----- SPES3. CARTE DE CHAMP CARTESIENNE MESUREE DU SPES3
C          D'APRES W. ROSCH, 1991
 47   CONTINUE
      KALC =2
      KUASEX = 3
      IF(READAT) CALL RCARTE(I1,I2,
     >                             ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- CARTE DE Champ CARTESIENNE 2-D DE CHALUT
 48   CONTINUE
      KALC =2
      KUASEX = 4
      IF(READAT) CALL RCARTE(I1,I2,
     >                             ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C-----  SPNTRK. Switch spin tracking
 49   CONTINUE
      IF(READAT) CALL RSPN
      IF(FITGET) CALL FITGT1
      CALL SPN
      GOTO 998
C----- SPNPRT. PRINT SPIN STATES
 50   CONTINUE
      CALL SPNPRT(LABEL(NOEL,1))
      GOTO 998
C----- BEND. Dipole magnet
 51   CONTINUE
      KALC =1
      KUASEX = 8
      IF(READAT) THEN 
        CALL RBEND(ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- SOLENOID
 52   CONTINUE
      KALC = 3
      KUASEX = 20
      IF(READAT) THEN 
        CALL RSOLEN(NDAT,NOEL,MXL,A,ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL  QUASEX(ND(NOEL))
      GOTO 998
C----- Plot transverse coordinates (for use on a workstation)
 53   CONTINUE
      GOTO 998
C----- SPNPRNL. Store  state of spins in logical unit
 54   CONTINUE
      IF(READAT) CALL RSPNST(I0,
     >                         PRLBSP,KPRTSP,LBLSP,NLB)
      IF(TA(NOEL,1).NE.'none') THEN
        NLB = 0
        TXTEMP = TA(NOEL,1)
        TXTEMP=TXTEMP(DEBSTR(TXTEMP):FINSTR(TXTEMP))
        CALL SPNPRW(TXTEMP,LBLSP,NLB)
        CALL SPNPRN(I0,NOEL-1,KLE(IQ(NOEL-1)),LABEL(NOEL-1,1),
     >                                             LABEL(NOEL-1,2))
      ENDIF
      GOTO 998
C----- CAVITE ACCELERATRICE
 55   CONTINUE
      IF(READAT) CALL RCAVIT
      IF(FITGET) CALL FITGT1
      CALL CAVITE(*999)
      GOTO 998
C----- PARTICUL.  DATA PARTICLE
 56   CONTINUE
C     .... MASSE(MeV/c2), CHARGE(C), G-yromagn., com life time(s), Free
      IF(READAT) CALL RPARTI(NDAT,NOEL,
     >                                 A)
      IF(FITGET) CALL FITGT1
      CALL PARTIC
      GOTO 998
C----- COMMANDE DES SCALINGS ( B, FREQ...)
 57   CONTINUE
      IF(READAT) CALL RSCAL
      IF(FITGET) CALL FITGT1
      CALL SCALIN
      GOTO 998
C----- ELREVOL.  ELCTROSTATIQ  FIELD  Ex(R=0,X) 1-D MEASURED.  CYLINDRICAL SYMM.
 68   CONTINUE
      KFLD=LC
C----- BREVOL. 1-D field B(r=0,x) on-axis field map, with cylindrical symmetry. 
 58   CONTINUE
      KALC =2
      KUASEX = 8
      IF(READAT) CALL RCARTE(I1,I1,
     >                             ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- POISSON. CARTE DE Champ CARTESIENNE 2-D FABRIQUEE PAR POISSON
 59   CONTINUE
      KALC =2
      KUASEX = 5
      IF(READAT) CALL RCARTE(I1,I2,
     >                             ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- DIPOLE. Similar to DIPOLES, but with a single magnet
 60   CONTINUE
      KALC =1
      KUASEX = 31
      IF(READAT) CALL RDIP(ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL AIMANT(ND(NOEL))
      GOTO 998
C----- CARTE MESUREE SPECTRO KAON GSI (DANFISICS)
 61   CONTINUE
      KALC =2
      KUASEX = 6
      IF(READAT) CALL RCARTE(I1,I2,
     >                             ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- MAP2D. 2D B-FIELD MAP, NO SPECIAL SYMMETRY (P. Akishin, 07/1992)
 62   CONTINUE
      IF(ID.LT.3)
     >   CALL ENDJOB('Use of MAP2D :  you need ID=3 in PARIZ.H',-99)
      KALC =2
      KUASEX = 9
      IF(READAT) CALL RCARTE(I1,I2,
     >                             ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- B DECAPOLE ET DERIVEES CALCULES EN TOUT POINT (X,Y,Z)
 63   CONTINUE
      KALC = 3
      KUASEX = 5
      IF(READAT) THEN
        CALL RQSOD(NDAT,NOEL,MXL,A,ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- B DODECAPOLE ET DERIVEES CALCULES EN TOUT POINT (X,Y,Z)
 64   CONTINUE
      KALC = 3
      KUASEX = 6
      IF(READAT) THEN
        CALL RQSOD(NDAT,NOEL,MXL,A,ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- WIENFILT - INTEGR NUMERIQ.
 65   CONTINUE
      KALC = 3
      KUASEX = 21
      KFLD=ML
      IF(READAT) CALL RWIENF(ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- FAISTORE - Similar to FAISCNL, with additional options:
C        - Print after LABEL'ed elements
C        - Print every other IPASS = mutltiple of IA
 66   CONTINUE
      IF(READAT) CALL RFAIST(MLB,
     >                           PRLB,KPRT,LBL,NLB)
      TXTEMP = TA(NOEL,1)
      TXTEMP=TXTEMP(DEBSTR(TXTEMP):FINSTR(TXTEMP))
      IF(TXTEMP.NE.'none') THEN
        IF(TA(NOEL,2).NE.'none') THEN
          CALL IMPFAW(TXTEMP,LBL,NLB)
          IF(.NOT. PRLB) CALL IMPFAI(KPRT,NOEL-1,KLE(IQ(NOEL-1)),
     >                           LABEL(NOEL-1,1), LABEL(NOEL-1,2))
        ENDIF
      ENDIF
      GOTO 998
C----- SPNSTORE - Similar to SPNPRNL, with additional options:
C        - Print after LABEL'ed elements
C        - Print every other IPASS = mutltiple of IA
 67   CONTINUE
      IF(READAT) CALL RSPNST(MLB,
     >                           PRLBSP,KPRTSP,LBLSP,NLB)
      TXTEMP = TA(NOEL,1)
      TXTEMP=TXTEMP(DEBSTR(TXTEMP):FINSTR(TXTEMP))
      IF(TXTEMP.NE.'none') THEN
        IF(TA(NOEL,2).NE.'none') THEN
          CALL SPNPRW(TXTEMP,LBLSP,NLB)
          IF(.NOT. PRLBSP) CALL SPNPRN(KPRTSP,NOEL-1,KLE(IQ(NOEL-1)),
     >                           LABEL(NOEL-1,1), LABEL(NOEL-1,2))
        ENDIF
      ENDIF
      GOTO 998
C----- EL2TUB - LENTILLE ELECTROSTATIQ A 2 TUBES OU DIAPHRAG.
 69   CONTINUE
      KALC = 3
      KUASEX = 22
      KFLD=LC
      IF(READAT) CALL REL2TU(ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- UNIPOT. LENTILLE ELECTROSTATIQ A 3 TUBES OU DIAPHRAG.
 70   CONTINUE
      KALC = 3
      KUASEX = 23
      KFLD=LC
      IF(READAT) CALL RUNIPO(ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- ELMULT. E MULTIPOLAIRE ET DERIVEES CALCULES EN TOUT POINT (X,Y,Z)
 71   CONTINUE
      KALC = 3
      KUASEX = MPOL+1
      KFLD=LC
      IF(READAT) THEN
        CALL RMULTI(NDAT,NOEL,MXL,A,MPOL,ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- EBMULT. E+B MULTIPOLAIRES ET DERIVEES CALCULES EN TOUT POINT (X,Y,Z)
 72   CONTINUE
      KALC = 3
      KUASEX = MPOL+1
      KFLD=ML
      IF(READAT) CALL REBMUL(ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- ESPACE LIBRE VIRTUEL
 73   CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1),A(NOEL,2)
      IF(FITGET) CALL FITGT1
      CALL ESLVIR(1,ZERO,ZERO)
      GOTO 998
C----- TRANSFORMATION FO(I,n)= A(,3)*FO(I,n) + A(,4)*FO(J,n) + A(,5)
 74   CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1),A(NOEL,2),
     >           A(NOEL,3),A(NOEL,4),A(NOEL,5)
      IF(FITGET) CALL FITGT1
      CALL TROBJ(1)
      GOTO 998
C----- TRANSLATION  F(2,n)=F(2,n)+A(,1) et F(4,n)=F(4,n)+A(,2)
 75   CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1),A(NOEL,2)
      IF(FITGET) CALL FITGT1
      CALL TRAOBJ(1)
      GOTO 998
C----- POLARMES. CARTE DE Champ POLAIRE. Field map, cylindrical coordinates. 
 76   CONTINUE
      KALC =2
      KUASEX = 22
      IF(READAT) CALL RCARTE(I2,I2,
     >                             ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL AIMANT(ND(NOEL))
      GOTO 998
C----- TRANSLATION X,Y,Z ET ROTATION RX,RY,RZ
 77   CONTINUE
      IF(READAT) READ(NDAT,*) (A(NOEL,II),II=1,6)
      IF(FITGET) CALL FITGT1
      CALL TRAROT(A(NOEL,1),A(NOEL,2),A(NOEL,3)
     >  ,A(NOEL,4),A(NOEL,5),A(NOEL,6))
      GOTO 998
C----- SRLOSS. Switch synchrotron ligth on
 78   CONTINUE
      IF(READAT) CALL RSRLOS
      IF(FITGET) CALL FITGT1
      CALL SRLOSS(IMAX,IPASS,*999)
      GOTO 998
C----- PICKUPS. 
 79   CONTINUE
      IF(READAT) CALL RPCKUP
      IF(FITGET) CALL FITGT1
      CALL PICKUP
      GOTO 998
C----- OPTICS. Transport the beam matrix and print/store it after keyword[s].
 80   CONTINUE
      IF(READAT) THEN
        READ(NDAT,fmt='(A)') TXTEMP
        READ(TXTEMP(DEBSTR(TXTEMP):FINSTR(TXTEMP)),
     >  *,ERR=801,END=801) KOPTCS, LBLOPT
        READ(TXTEMP(DEBSTR(TXTEMP):FINSTR(TXTEMP)),
     >  *,ERR=803,END=803) KOPTCS, TXT1, KOPIMP
      ENDIF
      GOTO 802
 801  CONTINUE
      LBLOPT = 'all'
      GOTO 802
 803  CONTINUE
      KOPIMP = 0
      IF(NRES.GT.0) THEN
        WRITE(NRES,*)  ' '
        WRITE(NRES,*)  ' OPTICS keyword. EOF or ERR while reading'
     >  ,'KOPTCS, LBLOPT, KOPIMP from  zgoubi.dat' 
      ENDIF
 802  CONTINUE
      IF (KOPTCS .NE. 1) KOPTCS = 0
      IF(NRES.GT.0) THEN
        WRITE(NRES,*)  ' '
        WRITE(NRES,*)  ' KOPTCS =',KOPTCS,'  (off/on = 0/1)'
        WRITE(NRES,*)  ' '
        WRITE(NRES,*)  ' LBLOPT = ' //
     >  LBLOPT(DEBSTR(LBLOPT):FINSTR(LBLOPT))  //
     >  '  (will print beam matrix into zgoubi.res at those labels)' 
        WRITE(NRES,*)  ' '
        WRITE(NRES,*)  ' KOPIMP=',KOPIMP,
     >  '  (print betas into zgoubi.OPTICS.out, no/yes = 0/1)'
        WRITE(NRES,*)  ' '
      ENDIF
C      OKLNO = .FALSE. 
      IF(KOPIMP.EQ.1) THEN
        IF(.NOT. OKLNO) THEN
          IF(IDLUNI(
     >          LNOPTI)) THEN
            OPEN(UNIT=LNOPTI,FILE='zgoubi.OPTICS.out',ERR=899)
            OKLNO = .TRUE.
          ENDIF
          IF(OKLNO) THEN
            WRITE(LNOPTI,FMT='(A)') '# FROM OPTICS KEYWORD'
            WRITE(LNOPTI,FMT='(A)') 
     >         '# ALFX,         BTX,          ALFY,         BTY, ' //
     >         '         ALFL,         BTL,          DX,         ' //
     >         '  DXP,          DY,           DYP,          PHIX/' //
     >         '2PI,     PHIY/2PI,     SUM_S,        #LMNT, X,   ' //
     >         '         XP,           Y,            YP,         ' //
     >         'KEYWORD,   LABEL1,    LABEL2       FO(6,1)       ' //
     >         'K0*L          K1*L          K2*L  '
            WRITE(LNOPTI,FMT='(A)') 
     >         '# 1             2             3             4    ' //
     >         '         5             6             7           ' //
     >         '  8             9             10            11   ' //
     >         '         12            13            14     15   ' //
     >         '         16            17            18          ' //
     >         '19         20         21           22            ' //
     >         '23            24            25      '
          endif
        endif
      ELSE
        OKLNO = .FALSE.
      ENDIF
      CALL OPTIC2(
     >            OKLNO,LNOPTI)
      GOTO 998
 899  CONTINUE
C      IF(NRES.GT.0) 
      WRITE(abs(NRES),*) 'KEYWORD OPTICS. Error open zgoubi.OPTICS.out.'
      GOTO 998
C-----  GASCAT. Switch gas-scattering
 81   CONTINUE
      IF(READAT) CALL RGASCA
      IF(FITGET) CALL FITGT1
      CALL GASINI(*999)
      GOTO 998
C----- UNDULATO.  UNDULATOR
 82   CONTINUE
      KALC =3
      KUASEX = 30
      IF(READAT) CALL RUNDUL(ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- ELCYLDEF
 83   CONTINUE
      IF(READAT) CALL RELCYL(ND(NOEL))
      IF(FITGET) CALL FITGT1
      KALC = 3
      KUASEX = 24
      KFLD=LC
      CALL AIMANT(ND(NOEL))
      GOTO 998
C----- ELMIR
 84   CONTINUE
      KALC = 3
      KUASEX = 25
      KFLD=LC
      IF(READAT) THEN
        CALL RELMIR(ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- ELCMIR
 85   CONTINUE
      KALC = 3
      KUASEX = 26
      KFLD=LC
      IF(READAT) THEN
        CALL RELCMI(ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL AIMANT(ND(NOEL))
      GOTO 998
C----- MAP2D-E: 2D E-FIELD MAP
 86   CONTINUE
      KFLD=LC
      GOTO 62
C----- SRPRNT. Print/Store S.R. loss tracking statistics into file
 87   CONTINUE
      CALL SRPRN(I0,NRES,IMAX)
      GOTO 998
C----- BETATRON. Betatron core
 88   CONTINUE
      IF(READAT) READ(NDAT,*) DPKCK
      IF(FITGET) CALL FITGT1
      CALL DPKICK(DPKCK)
      GOTO 998
C----- TWISS. Compute linear lattice functions, chromaticities, etc. 
C      Also prints periodic beta functions (by setting KOPTCS to 1).
 89   CONTINUE
C                            ktwiss=1 :  Fac_dp   Fac-ampl
C                            ktwiss=2 :  Prtcl#   unused    [coupled]
      NOELS = NOEL
      IF(READAT) THEN
        READ(NDAT,fmt='(a)') TXT132
        IF(STRCON(TXT132,'!',
     >                       IS)) TXT132 = TXT132(1:IS-1)
        IF(STRCON(TXT132,'coupled',
     >                             IS)) THEN
          OKCPLD = .TRUE. 
          TA(NOEL,1) = 'coupled'
        ENDIF
        READ(TXT132,*) A(NOEL,1),A(NOEL,2),A(NOEL,3)
      ENDIF
C      OKCPLD = 
C     >TA(NOEL,1)(DEBSTR(TA(NOEL,1)):FINSTR(TA(NOEL,1))) .eq. 'coupled'
      IF(.NOT. OKLNO) THEN
        IF(IDLUNI(
     >            LNOPTI)) THEN
          OPEN(UNIT=LNOPTI,FILE='zgoubi.TWISS.out',ERR=899)
          OKLNO = .TRUE.
        ENDIF
      ENDIF
      CALL OPTIC2(
     >            OKLNO,LNOPTI)

      CALL TWISS(LNOPTI,OKCPLD,
     >                  KOPTCS,READAT,KTW)
      IF(KOPTCS .EQ. 1) THEN
        KOPIMP = 2
        LBLOPT = 'all'
        PRDIC = .TRUE.
          IF(OKLNO) THEN
            WRITE(LNOPTI,fmt='(a)') '# From TWISS keyword'
            WRITE(LNOPTI,fmt='(a)') 
     >         '# alfx,         btx,          alfy,         bty, ' //
     >         '         alfl,         btl,          Dx,         ' //
     >         '  Dxp,          Dy,           Dyp,          phix/' //
     >         '2pi,     phiy/2pi,     sum_s,        #lmnt, x,   ' //
     >         '         xp,           y,            yp,         ' //
     >         'KEYWORD,   label1,    label2       FO(6,1)       ' //
     >         'K0*L          K1*L          K2*L          |C|    '
            WRITE(LNOPTI,fmt='(a)') 
     >         '# 1             2             3             4    ' //
     >         '         5             6             7           ' //
     >         '  8             9             10            11   ' //
     >         '         12            13            14     15   ' //
     >         '         16            17            18          ' //
     >         '19         20         21           22            ' //
     >         '23            24            25            26 '
          ENDIF
      ELSE
C        KOPIMP = 0
C        OKLNO = .FALSE.
      ENDIF
      IF(IPASS .EQ. 2*KTW+1) READAT = .TRUE.
c      write(*,*) ' zgoubi, noel, ktw, readat : ',noel,ktw,ipass,readat
c      write(*,*) ' zgoubi, noel, ktw, readat : ',noel,ktw,ipass,readat
      GOTO 998
C----- END. End of run, except for some options that may need more
 90   CONTINUE
      CALL END(
     >         READAT,NOEL,*998)
      GOTO 999
C----- FFAG. FFAG Sector multi-dipole. 
 91   CONTINUE
      KALC = 1
      KUASEX = 27
      IF(READAT) THEN
        CALL RFFAG(ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL AIMANT(ND(NOEL))
      GOTO 998
C----- HELIX. Helical field (twisted dipole)
 92   CONTINUE
      KALC = 3
      KUASEX = 28
      IF(READAT) THEN
        CALL RHELIX(ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C-----  CSR. Switch  coherent SR interaction
 93   CONTINUE
      IF(READAT) CALL RCSR
      IF(FITGET) CALL FITGT1
      CALL CSRI
      GOTO 998
C----- Each particle is pushed a Delta-S distance
 94   CONTINUE
      IF(READAT) READ(NDAT,*) A(NOEL,1)
      IF(FITGET) CALL FITGT1
      CALL PATH(A(NOEL,1))
      GOTO 998
C----- COILS. 
 95   CONTINUE
      KALC = 3
      KUASEX = 29
      IF(READAT) THEN 
        CALL RCOILS(NDAT,NOEL,MXL,A,ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL  QUASEX(ND(NOEL))
      GOTO 998
C----- GETFITVAL.  Get parameter values resulting from FIT, stored in TA(NOEL,1)
 96   CONTINUE
        CALL FITSTA(I5,
     >                 FITING)
C      IF(READAT) THEN
      IF(.NOT.FITING) THEN
        if(ipass.eq.1)CALL RFITGT
        CALL FITGTV(TA(NOEL,1),
     >                         FITGET)
      ENDIF
C          WRITE(*,*) ' ZGOUBI   GETFITVAL ',FITING,NOEL,IPASS
      GOTO 998
C----- SUPERPOSE.  To superimpose magnets. 2B developped
 97   CONTINUE
C      IF(READAT) CALL RSUPER(NDAT,NOEL,MXL,A)
      KUASEX = 31
      IF(FITGET) CALL FITGT1
C      CALL SUPERP
      GOTO 998
C----- MARKER. 
 98   CONTINUE
      KUASEX = 0  
      IF(LABEL(NOEL,2).EQ.'.plt' .OR. LABEL(NOEL,2).EQ.'.PLT') THEN
        CALL OPEN2('MAIN',NPLT,FILPLT)
        DO IT = 1, IMAX
          TRAD = F(3,IT)*1.D-3
          PRAD = F(5,IT)*1.D-3
          CALL IMPPLB(NPLT,F(2,IT),TRAD,F(4,IT),PRAD,ZERO,
     >    F(6,IT),F(7,IT),ZERO,AMQ(1,IT),AMQ(2,IT),IEX(IT),IT)
        ENDDO
      ENDIF
      GOTO 998
C----- DIPOLES. A set of neiboring or overlapping dipoles. 
 99   CONTINUE
      KALC =1
      KUASEX = 32
      IF(READAT) CALL RDIPS(ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL AIMANT(ND(NOEL))
      GOTO 998
C----- TRACKING. 
 100  CONTINUE
        READ(NDAT,*) NLMA, NLMB 
        WRITE(NRES,*) ' Tracking,  NLM_A  ->  NLM_B : ',NLMA,' -> ',NLMB
        CALL TRACK(NLMA,NLMB)
        CALL ENDJOB(' End of job after TRACKING',-99)
      GOTO 998
C----- FFAG-SPI. FFAG, spiral. 
 101  CONTINUE
      KALC = 1
      KUASEX = 33
      IF(READAT) THEN
        CALL RFFAG(ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL AIMANT(ND(NOEL))
      GOTO 998
C----- EMMA. Read  2-D or 3-D field map, 
C      with mesh either cartesian (KART=1) or cylindrical (KART=2). 
 103  CONTINUE
      KALC = 2
      IF(READAT) CALL REMMA(KART,I3,
     >                              ND(NOEL))
      IF    (A(NOEL,22) .EQ. 1) THEN
C------- 2-D map, KZMA = 1
        KUASEX = 34
      ELSEIF(A(NOEL,22) .GT. 1) THEN
C------- 3-D map, KZMA > 1
        KUASEX = 35
        IF(IZ.LE.1) CALL ENDJOB(' *** ERROR ; cannot use a 3-D map, need
     >  recompile zgoubi, using IZ>1 in PARIZ.H',-99) 
      ENDIF
      IF(FITGET) CALL FITGT1
      IF    (KART.EQ.1) THEN 
        CALL QUASEX(ND(NOEL))
      ELSEIF(KART.EQ.2) THEN 
        CALL AIMANT(ND(NOEL))
      ENDIF
      GOTO 998
C----- DIPOLEC. Like DIPOLES, with cartesian coordinates
 104  CONTINUE
      KALC =1
      KUASEX = 36
      IF(READAT) CALL RDIPC(ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- REVERSE. 
 105  CONTINUE
      CALL REVERS
      GOTO 998
C----- SYSTEM. System call
 106  CONTINUE
C        stop '***************'
C      IF(READAT) READ(NDAT,*) A(NOEL,1)
      READ(NDAT,*) A(NOEL,1)
      NCMD = NINT(A(NOEL,1))
      IF(NRES.GT.0) WRITE(NRES,*) ' Number of commands : ',NCMD
C      IF(READAT) THEN 
        DO I = 1, NCMD
          READ(NDAT,FMT='(A)') SYSCMD
          CALL SYSTEM(SYSCMD(DEBSTR(SYSCMD):FINSTR(SYSCMD)))
          IF(NRES.GT.0) WRITE(NRES,*) 
     >    SYSCMD(DEBSTR(SYSCMD):FINSTR(SYSCMD))
        ENDDO 
C      ENDIF
          WRITE(*,*) 
     >    SYSCMD(DEBSTR(SYSCMD):FINSTR(SYSCMD))
                  read(*,*)
      GOTO 998
C----- SPINR. Spin rotator 
 107  CONTINUE
      IF(READAT) CALL RSPINR
      IF(FITGET) CALL FITGT1
      CALL SPINR
      GOTO 998
C----- BENDTH. Pure dipole field, analytical push. 
 108  CONTINUE
      IF(READAT) CALL RBNDTH(
     >                       ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL BNDTHI(ND(NOEL))
      GOTO 998
C----- AGSMM. AGS main magnet. Works like MULTIPOL + various refinements or specificities. 
 109  CONTINUE
      KALC = 3
      KUASEX = 37
      IF(READAT) THEN
        CALL RAGSMM(NDAT,NOEL,MXL,A,ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- BEAMBEAM.
 110  CONTINUE
      IF(READAT) CALL RBB
      IF(FITGET) CALL FITGT1
      CALL BB
      GOTO 998
C----- AGSQUAD. AGS quadrupole. It has 2 windings, works otherwise like QUADRUPO.
 111  CONTINUE
      KALC = 3
      KUASEX = 38
      IF(READAT) THEN
        CALL RAGSQU(NDAT,NOEL,MXL,A,ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- SYNRAD.
 112  CONTINUE
      CALL SYNTRK
      GOTO 998
C----- OPTION.
 113  CONTINUE
      IF(READAT) CALL ROPTIO
      CALL OPTION
      GOTO 998
C----- EPLATES. 
 114  CONTINUE
      KALC = 3
      KUASEX = 39
      IF(READAT) THEN 
        CALL REPLAT(ND(NOEL))
      ELSE
        CALL STPSI1(NOEL)
      ENDIF
      IF(FITGET) CALL FITGT1
      CALL QUASEX(ND(NOEL))
      GOTO 998
C----- DAMPER. 
 115  CONTINUE
      IF(READAT) CALL RDAMPE
      IF(FITGET) CALL FITGT1
      CALL DAMPER
      GOTO 998
C----- CYCLOTRON.
 116  CONTINUE
      KALC =1
      KUASEX = 40
      IF(READAT) CALL RCYCLO(ND(NOEL))
      IF(FITGET) CALL FITGT1
      CALL AIMANT(ND(NOEL))
      GOTO 998
C----- ERRORS. 
 117  CONTINUE
      IF(READAT) CALL RERROR
      IF(FITGET) CALL FITGT1
      CALL ERRORS
      GOTO 998

C-------------------------
C-------------------------
      ENTRY ZGLMNT(
     >             TXTELO)
      TXTELO = TXTELT
      RETURN
C Current KLEY
      ENTRY ZGKLEY( 
     >             KLEO)
      KLEO = KLEY
      RETURN
C Size of KLE array
      ENTRY ZGMXKL( 
     >             MXKLEO)
      MXKLEO = MXKLE
      RETURN
      ENTRY ZGNOEL(
     >             NOELO)
C Current elmnt #
      NOELO = NOEL
      RETURN
C KLEY[IKL]
      ENTRY ZGKLE(IKL, 
     >                KLEO)
      KLEO = KLE(IKL)
      RETURN
      ENTRY ZGMXL( 
     >            NBLMNO)
      NBLMNO = NBLMN
      RETURN
      ENTRY ZGPNLT( 
     >             PNLTGO)
      PNLTGO = PNLTGT
      RETURN
      ENTRY ZGIPAS( 
     >             IPASSO)
C Current pass #
      IPASSO = IPASS
      RETURN

      END
