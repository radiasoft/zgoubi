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
      SUBROUTINE DIPSI(SCAL,
     >                      DSREF,IRD,IDB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.AIM_2.H"     ! COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.DROITE_2.H"     ! COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
C      INCLUDE "C.ORDRES.H"     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI

      DIMENSION FTAB(5,5)

      PARAMETER (NM=5)
      DIMENSION NBFACE(NM)
      DIMENSION CE(NM,6),CS(NM,6),CC(NM,6)
      DIMENSION ACN(NM),DRM(NM),RRM(NM),HNORM(NM)
      DIMENSION IND(NM), CIND(NM,10)
      DOUBLE PRECISION LAMBDE, LAMBDS, LAMBD3
      DIMENSION LAMBDE(NM),QSIE(NM),NCOEFE(NM),SHIFTE(NM)
      DIMENSION LAMBDS(NM),QSIS(NM),NCOEFS(NM),SHIFTS(NM)
      DIMENSION LAMBD3(NM),QSI3(NM),NCOEF3(NM),SHIFT3(NM),RM3(NM)
      DIMENSION UMEGA(NM),THETA(NM),R1(NM),U1(NM),U2(NM),R2(NM)
      DIMENSION UMEGAS(NM),THETAS(NM),R1S(NM),U1S(NM),U2S(NM),R2S(NM)
      DIMENSION UMEGA3(NM),THETA3(NM),R13(NM),U13(NM),U23(NM),R23(NM)

      SAVE ACN,DRM,HNORM,CE,CS,CC
      SAVE IND,CIND
      SAVE LAMBDE,QSIE,NCOEFE,SHIFTE
      SAVE LAMBDS,QSIS,NCOEFS,SHIFTS
      SAVE LAMBD3,QSI3,NCOEF3,SHIFT3
      SAVE UMEGA,THETA,R1,U1,U2,R2
      SAVE UMEGAS,THETAS,R1S,U1S,U2S,R2S
      SAVE NN, RESOL, ITYPF
      SAVE NBMAG, NBFACE

      DIMENSION BZ0(5,5)

      CHARACTER(14) TYPCAL(2)

      SAVE TYPCAL

      PARAMETER (PLIM=80.D0)
      LOGICAL SHARPE, SHARPS

      SAVE FINTE, FINTS

      SAVE AD20, AD22, AD24, AD26, AD28
     >   , AD30, AD32, AD34, AD36, AD38
     >   , AD40, AD42, AD44, AD46, AD48
     >   , AD50, AD52, AD54, AD56, AD58
      SAVE AF20, AF22, AF24, AF26, AF28
     >   , AF30, AF32, AF34, AF36, AF38
     >   , AF40, AF42, AF44, AF46, AF48
     >   , AF50, AF52, AF54, AF56, AF58

      DATA TYPCAL / ' analytic', ' interpolation'/

C  Isochronous e-model, Rees
C      DATA FINTE, FINTS / 0.025D0, 0.025D0 /
C  Isochronous muon FFAG, Rees
C      DATA FINTE, FINTS / 0.4D0, 0.4D0 /
      DATA FINTE, FINTS / 0.2D0, 0.2D0 /

C Sandro's ffag field
      DATA AD20, AD22, AD24, AD26, AD28 /
     >597.5663294408266D0    , -9.380049855585632D0,
     >0.02312190362096569D0  ,-0.000016305002117346734D0,
     >3.4992136271753898D-9/
      DATA AD30, AD32, AD34, AD36, AD38 /
     >-1.4342562745299526D6, 12959.155925600144D0,
     >-38.24911496883343D0   , 0.037112611054756564D0,
     >-0.000011241547359038D0/
      DATA AD40, AD42, AD44, AD46, AD48 /
     >2.3601804363826987D8 , -3.288355034691052D6,
     > 10335.690227105812D0  , -10.536007335638033D0,
     > 0.003307235670320862D0/
      DATA AD50, AD52, AD54, AD56, AD58 /
     >-2.022113627491819D10, 2.718621791696837D8,
     >-859226.1506962734D0   , 893.6985655822829D0,
     >-0.28514145252461576D0 /

      DATA AF20, AF22, AF24, AF26, AF28 /
     >784.862D0, 34.9238D0, -0.873452D0,  0.0056483D0, -0.0000105406D0/
      DATA AF30, AF32, AF34, AF36, AF38 /
     >-20777.7D0,  8027.21D0, -177.955D0, 1.12474D0,-0.00208022D0/
      DATA AF40, AF42, AF44, AF46, AF48 /
     >-1.96463034586771D6, 462200.886608188D0, -9900.4400446493D0,
     >62.1655878585100D0,-0.11469175008766148D0/
      DATA AF50, AF52, AF54, AF56, AF58 /
     >-7.99910758615139D7, 9.16588708380287D6,
     >-186561.6360278200D0, 1156.849173782191D0,-2.1235545491928156D0/

C  NBMAG=number of magnets.  AT=total extent angle of field
      NP = 2
      NBMAG = NINT(A(NOEL,NP))
      NP=NP+1
      AT    = A(NOEL,NP)
      NP=NP+1
      RM    = A(NOEL,NP)

      IF(NRES.GT.0) WRITE(NRES,200) NBMAG,AT,RM
  200 FORMAT(20X,1P,'Dipole  magnet  N-tuple. ',
     >'Number  of  sectors  N = ',I2,//,
     >11X,' Total  angular  extent  of  sector  AT = ',G12.4,' deg.',
     >/,11X,' Reference  radius  of  sector  RM = ',G12.4,' cm')

C  HNORM=Champ MAX DANS LE DIPOLE.
C  CIND= field indices
C  AT=ANGLE TOTAL DE LA zone DE Champ, ACENT='ANGLE AU CENTRE',
C    RM=RAYON MOYEN DE LA zone DE Champ.
C  NBFACE(KMAG)=(2)3 : dipole limited by (2)3 field boundaries

C-----  a list of NBMAG.LE.NM magnets
      KMAG = 0
 10   CONTINUE
      KMAG = KMAG+1
      NP=NP+1
      ACENT = A(NOEL,NP)
      NP=NP+1
      DRM(KMAG)    = A(NOEL,NP)
      NP=NP+1
      HNORM(KMAG) = A(NOEL,NP)*SCAL
      NP=NP+1
      IND(KMAG) = NINT(A(NOEL,NP))
      IF(IND(KMAG).NE.0) THEN
        INDX = IND(KMAG)
        IF(INDX.LT.0) INDX = -INDX
        IF(INDX.GT.10) STOP ' *** SBR DIPSI, IND cannot be > 10'
        DO 19 KND = 1, INDX
          NP=NP+1
          CIND(KMAG,KND) = A(NOEL,NP)
 19     CONTINUE
      ENDIF

      ACN(KMAG) = ACENT * RAD
      IF(NRES.GT.0) THEN
        WRITE(NRES,100) KMAG,ACN(KMAG)/RAD,DRM(KMAG),HNORM(KMAG)
 100    FORMAT(/,5X,'Dipole  magnet # ',I1,/,1P,
     >  11X,' Positionning  angle ACENT : ',G12.4,' degrees',/,
     >  11X,' Positionning  wrt.  RM  : ',G12.4,' cm',/,
     5  11X,' HNORM =',G12.4,' KGauss')
C,5X,'COEF.N =',F9.5,5X,'COEF.B =',
C     6  F9.5,5X,'COEF.G=',F9.5)
        IF(IND(KMAG).NE.0) THEN
          WRITE(NRES,101) (KND,CIND(KMAG,KND),KND=1,INDX)
 101      FORMAT(11X,' Field  index  ',I2,' :   ',1P,G12.4,2X)
        ELSE
          WRITE(NRES,FMT='(11X,'' Field  indices :  NONE'')')
        ENDIF
      ENDIF

C Entrance EFB
      NP=NP+1
      LAMBDE(KMAG) = A(NOEL,NP)
      NP=NP+1
      QSIE(KMAG) = A(NOEL,NP)
      NP=NP+1
      NCOEFE(KMAG) = NINT(A(NOEL,NP))
      DO 227 I=1,6
        NP=NP+1
 227    CE(KMAG,I) = A(NOEL,NP)
      NP=NP+1
      SHIFTE(KMAG) = A(NOEL,NP)

      NP=NP+1
      UMEGA(KMAG) = A(NOEL,NP)
      NP=NP+1
      THETA(KMAG) = A(NOEL,NP)
      NP=NP+1
      R1(KMAG)    = A(NOEL,NP)
      NP=NP+1
      U1(KMAG)    = A(NOEL,NP)
      NP=NP+1
      U2(KMAG)    = A(NOEL,NP)
      NP=NP+1
      R2(KMAG)    = A(NOEL,NP)

      IF(NRES.GT.0) THEN
        WRITE(NRES,106) 'Entrance',LAMBDE(KMAG),QSIE(KMAG)
106     FORMAT (/,5X,A8,'  face',/,10X,
     >  'Reference gap for fringe field : g_0 =',
     >  1P,G12.4,' CM ,  gap-k =',G12.4)
        WRITE(NRES,127) NCOEFE(KMAG),(CE(KMAG,I),I=1,6),SHIFTE(KMAG)
127     FORMAT (10X,' COEFFICIENTS :',I3,6F10.5
     2  ,/,10X,' Shift  of  EFB  = ',G12.4,' CM',/)
        IF(LAMBDE(KMAG) .EQ. 0.D0)  WRITE(NRES,FMT='(/,
     >      ''  ***  Warning : sharp edge '',
     >      ''model entails vertical wedge focusing approximated with'',
     >        '' first order kick  ***'')')
        WRITE(NRES,103) UMEGA(KMAG),THETA(KMAG),R1(KMAG),U1(KMAG),
     >          U2(KMAG),R2(KMAG)
 103    FORMAT(10X,'OMEGA =',F7.2,' deg.',
     >   5X,'Wedge  angle  =',F7.2,' deg.',/,
     1   11X,'Radius 1 =',1P,G10.2,' cm',/ ,
     2   11X,'Straight  segment 1 =',G10.2,' cm',/ ,
     3   11X,'Straight  segment 2 =',G10.2,' cm',/ ,
     4   11X,'Radius 2 =',G10.2,' cm')
        IF(R1(KMAG)*R2(KMAG) .EQ. 0.D0) WRITE(NRES,123)
 123    FORMAT(10('*'),' ATTENTION |    R1 OU R2 = 0 |',10('*'))
      ENDIF

C Exit EFB
      NP=NP+1
      LAMBDS(KMAG) = A(NOEL,NP)
      NP=NP+1
      QSIS(KMAG)   = A(NOEL,NP)
      NP=NP+1
      NCOEFS(KMAG) = NINT(A(NOEL,NP))
      DO 228 I=1,6
        NP=NP+1
 228    CS(KMAG,I) = A(NOEL,NP)
      NP=NP+1
      SHIFTS(KMAG) = A(NOEL,NP)

      NP=NP+1
      UMEGAS(KMAG)  = A(NOEL,NP)
      NP=NP+1
      THETAS(KMAG) = A(NOEL,NP)
      NP=NP+1
      R1S(KMAG)    = A(NOEL,NP)
      NP=NP+1
      U1S(KMAG)    = A(NOEL,NP)
      NP=NP+1
      U2S(KMAG)    = A(NOEL,NP)
      NP=NP+1
      R2S(KMAG)    = A(NOEL,NP)
      IF(NRES.GT.0) THEN
        WRITE(NRES,106) 'Exit    ',LAMBDS(KMAG),QSIS(KMAG)
        WRITE(NRES,127) NCOEFS(KMAG),(CS(KMAG,I),I=1,6),SHIFTS(KMAG)
        IF(LAMBDS(KMAG) .EQ. 0.D0)  WRITE(NRES,FMT='(/,
     >      ''  ***  Warning : sharp edge '',
     >      ''model entails vertical wedge focusing approximated with'',
     >        '' first order kick  ***'')')
        WRITE(NRES,103) UMEGAS(KMAG),THETAS(KMAG),R1S(KMAG),U1S(KMAG),
     >              U2S(KMAG),R2S(KMAG)
        IF(R1S(KMAG)*R2S(KMAG) .EQ. 0.D0) WRITE(NRES,123)
      ENDIF

C Lateral EFB
      NP=NP+1
      LAMBD3(KMAG) = A(NOEL,NP)
      NP=NP+1
      QSI3(KMAG)   = A(NOEL,NP)
      NP=NP+1
      NCOEF3(KMAG) = NINT(A(NOEL,NP))
      DO 229 I=1,6
        NP=NP+1
 229    CC(KMAG,I) = A(NOEL,49+I)
      NP=NP+1
      SHIFT3(KMAG) = A(NOEL,NP)

      NP=NP+1
      UMEGA3(KMAG) = A(NOEL,NP)
      NP=NP+1
      THETA3(KMAG) = A(NOEL,NP)
      NP=NP+1
      R13(KMAG)    = A(NOEL,NP)
      NP=NP+1
      U13(KMAG)    = A(NOEL,NP)
      NP=NP+1
      U23(KMAG)    = A(NOEL,NP)
      NP=NP+1
      R23(KMAG)    = A(NOEL,NP)
      NP=NP+1
      RM3(KMAG)    = A(NOEL,NP)

      NBFACE(KMAG)=2
      IF(QSI3(KMAG) .NE. 0.D0) NBFACE(KMAG)=3

      IF(NRES.GT.0) THEN
        IF(NBFACE(KMAG).EQ.3) THEN
          STOP ' Lateral EFB is not implemented. Use Entrance/Exit only'
          WRITE(NRES,106) 'Lateral ',LAMBD3(KMAG),QSI3(KMAG)
          WRITE(NRES,127) NCOEF3(KMAG),(CC(KMAG,I),I=1,6),SHIFT3(KMAG)
          WRITE(NRES,113) RM3(KMAG)
113       FORMAT(20X,' Face centred on direction ACENT+OMEGA,  at ',
     >    G12.4,' CM'/)
          WRITE(NRES,103) UMEGA3(KMAG),THETA3(KMAG),R13(KMAG),U13(KMAG),
     >         U23(KMAG),R23(KMAG)
          IF(R13(KMAG)*R23(KMAG) .EQ. 0.D0) WRITE(NRES,123)
        ENDIF
      ENDIF

      IF(KMAG.LT.NBMAG) GOTO 10
C-----------------------------

      SHARPE = .TRUE.
      SHARPS = .TRUE.
      DO 57 KMAG = 1, NBMAG
        SHARPE = SHARPE .AND. (LAMBDE(KMAG).EQ.0.D0)
        SHARPS = SHARPS .AND. (LAMBDS(KMAG).EQ.0.D0)
 57   CONTINUE
      IF(NRES.GT.0) WRITE(NRES,*) ' WARNING, SBR DIPSI : ',
     >' Make sure you want hard edge correction with ',
     >' FINT at entrance / exit = ', FINTE,FINTS
      IF(SHARPE) THEN
        GAPE = -LAMBDE(1)
        IF(NRES.GT.0)
     >  WRITE(NRES,FMT='(''Entrance hard edge is to be implemented'')')
        CALL INTEG1(ZERO,FINTE,GAPE)
      ENDIF
      IF(SHARPS) THEN
        GAPS = -LAMBDS(1)
        IF(NRES.GT.0)
     >  WRITE(NRES,FMT='(''Exit hard edge is to be implemented'')')
        CALL INTEG2(ZERO,FINTS,GAPS)
      ENDIF

C Get type of field & deriv. calculation
      NP=NP+1
      KIRD = NINT(A(NOEL,NP))
C Get resol, or idb
      NP=NP+1
      RESOL=A(NOEL,NP)
      IF    (KIRD.NE.0) THEN
C        interpolation method
        IF(SHARPE .OR. SHARPS) CALL ENDJOB
     >    ('ERROR :  sharp edge not compatible with num. deriv.',-99)
        IRD = KIRD
        KIRD = 1
        IF    (IRD.EQ.2) THEN
          NN=3
        ELSEIF(IRD.EQ.25) THEN
          NN=5
        ELSEIF(IRD.EQ.4) THEN
          NN=5
        ELSE
          STOP ' *** ERROR - SBR DIP2, WRONG VALUE IRD'
        ENDIF
      ELSEIF(KIRD.EQ.0) THEN
C        analytic
C        IDB = NINT(10*A(NOEL,NP))
        IDB = NINT(RESOL)
        IF(IDB.NE.4) IDB=2
      ENDIF
      CALL CHAMC6(KIRD)

C Get flag field type. iord.option has the form xx.yy, iord=2, 25 or 4 and yy=01-99
C      ITYPF = NINT(100.D0* (A(NOEL,NP-1)-DBLE(IRD)) )
      ITYPF = NINT(100.D0* (A(NOEL,NP-1)-DBLE(NINT(A(NOEL,NP-1)))) )
      IF(NRES .GT. 0) THEN
        WRITE(NRES,FMT='(/,5X,'' itypf= '',i6)') itypf
        IF    (ITYPF.EQ.0) THEN
          WRITE(NRES,FMT='(/,5X,'' Field formula : default  '')')
        ELSEIF(ITYPF.EQ.1) THEN
          WRITE(NRES,FMT='(/,5X,'' Field formula : type 2  '')')
        ELSEIF(ITYPF.EQ.7) THEN
C          Sandro's FFAG magnet
          WRITE(NRES,FMT='(/,5X,'' Field  coefficients, AD : '')')
          WRITE(NRES,FMT='(1P,5G14.6)')  AD20, AD22, AD24, AD26, AD28
     >    , AD30, AD32, AD34, AD36, AD38
     >    , AD40, AD42, AD44, AD46, AD48
     >    , AD50, AD52, AD54, AD56, AD58

        ELSEIF(ITYPF.EQ.8 .OR. ITYPF.EQ.9) THEN
          WRITE(NRES,FMT='(/,5X,'' Field  coefficients, AF : '')')
          WRITE(NRES,FMT='(1P,5G14.6)')  AF20, AF22, AF24, AF26, AF28
     >    , AF30, AF32, AF34, AF36, AF38
     >    , AF40, AF42, AF44, AF46, AF48
     >    , AF50, AF52, AF54, AF56, AF58
        ENDIF
      ENDIF

      AE=0.D0
      AS=0.D0
      AT = AT * RAD
      XI = 0.D0
      XF = AT

C Formula to be revisited...
      KMAG = 1
      DSREF = RM * (  (UMEGA(1)-UMEGAS(KMAG))*RAD +
     >         TAN(ACN(1) - UMEGA(1)* RAD) +
     >         TAN(AT - ACN(KMAG) + UMEGAS(KMAG)* RAD) )

      IF(NRES.GT.0) THEN
        WRITE(NRES,FMT='(/,5X,'' Field & deriv. calculation :'',A)')
     >  TYPCAL(KIRD+1)
        IF    (KIRD.NE.0) THEN
          IF(IRD .EQ. 2) THEN
            WRITE(NRES,121) '3*3', RESOL
 121        FORMAT(20X,A3,
     >        '-point  interpolation, size of flying mesh :  ',
     >        'integration step /',1P,G12.4)
          ELSEIF(IRD .GE. 4) THEN
C----------- IRD= 4 OR 25
            WRITE(NRES,121) '5*5', RESOL
          ENDIF
        ELSEIF(KIRD.EQ.0) THEN
          WRITE(NRES,FMT='(20X,
     >    ''Derivatives computed to order '',I1)') IDB
        ENDIF
      ENDIF

      RETURN

C--------------------------------------------------------------------
C----- Numerical calculation of B & derivatives ---------------------
      ENTRY DIPSF(TTA,RO,
     >                   DTTA,DRO,FTAB)
      CALL INTEG5(
     >            STEP)

      CALL RAZ(FTAB,5*5)
      DRO = STEP/RESOL
      DTTA = DRO/RM

      KMAG = 0
 20   CONTINUE
      KMAG = KMAG+1

      RRM(KMAG) = RM +DRM(KMAG)
C           write(*,*) ' dipsi rm, drm :',RM ,DRM(KMAG)

C Entrance EFB
C     ** POUR LES BESOINS DE LA SIMPLE PRECISION :
      IF(R1(KMAG)*R1(KMAG) .GE. 1.D6) U1(KMAG) = -1.D6
      IF(R2(KMAG)*R2(KMAG) .GE. 1.D6) U2(KMAG) =  1.D6

      UMEG = UMEGA(KMAG) * RAD
      TETA = THETA(KMAG) * RAD
      UT = UMEG - TETA
      SINO = SIN(UT)
      COSO = COS(UT)
      TANO = SINO / COSO

C  CALCUL DES PARAMETRES DE FACE  :
C  AXE X DU REFERENTIEL  = PARALLELE  DIRCTN ACENT, SENS DES RAYONS CROISSANTS
C  AXE Y DU REFERENTIEL  = ORTHOGONAL DIRCTN ACENT, DIRIGE VERS FACE D'ENTREE
C  ORIGINE DU REFERENTIEL= A L'INTERSECTION DE RM ET ACENT

C  PROJECTIONS DU RAYON MOYEN SUR LES AXES X,Y
      XB = RRM(KMAG) * ( COS(UMEG) - 1.D0)
      YB = RRM(KMAG) * SIN(UMEG)
C  COORDONNEES DE L'EXTREMITE A DE LA PARTIE LINEAIRE DE LONG. U2
      XA = U2(KMAG) * COSO + XB
      YA = YB + U2(KMAG) * SINO
C  COORDONNEES DE L'EXTREMITE C DE LA PARTIE LINEAIRE DE LONG. U1
      XC = U1(KMAG) * COSO + XB
      YC = YB + U1(KMAG) * SINO
C  COORDONNEES DU CENTRE DE COURBURE DE RAYON R1
      XD = XC + R1(KMAG) * SINO
      YD = YC - R1(KMAG) * COSO
C  COORDONNEES DU CENTRE DE COURBURE DE RAYON R2
      XE = XA + R2(KMAG) * SINO
      YE = YA - R2(KMAG) * COSO
      SIN2 = SINO**2
      COS2 = COSO**2
      SICO = SINO * COSO

C Exit EFB
C     ** POUR LES BESOINS DE LA SIMPLE PRECISION :
      IF(R1S(KMAG)*R1S(KMAG) .GE. 1.D10) U1S(KMAG) = -1.D6
      IF(R2S(KMAG)*R2S(KMAG) .GE. 1.D10) U2S(KMAG) =  1.D6

      UMEGS = UMEGAS(KMAG) * RAD
      TETAS = THETAS(KMAG) * RAD
      UTS= UMEGS- TETAS
      SINOS= SIN(UTS)
      COSOS= COS(UTS)
      TANOS= SINOS/ COSOS
      XBS = RRM(KMAG) * ( COS(UMEGS) - 1.D0)
      YBS = RRM(KMAG) * SIN(UMEGS)
      XAS = U2S(KMAG) * COSOS + XBS
      YAS = YBS + U2S(KMAG) * SINOS
      XCS = U1S(KMAG) * COSOS + XBS
      YCS = YBS + U1S(KMAG) * SINOS
      XDS = XCS + R1S(KMAG) * SINOS
      YDS = YCS - R1S(KMAG) * COSOS
      XES = XAS + R2S(KMAG) * SINOS
      YES = YAS - R2S(KMAG) * COSOS
      SIN2S = SINOS**2
      COS2S = COSOS**2
      SICOS = SINOS * COSOS

C Third EFB
      IF(NBFACE(KMAG) .EQ. 3) THEN
        UMEG3 = UMEGA3(KMAG) * RAD
        TETA3 = THETA3(KMAG) * RAD
        UT3= UMEG3- TETA3
        IF(ABS(UT3) .GE. .5d0*PI) THEN
          IF(NRES.GT.0) WRITE(NRES,139)
139       FORMAT(/,20X,' EFB lateral, ATTENTION : ABS(UT3) > 90 deg =>',
     >    '  ambiguous  for  calculation  of  field')
        ENDIF

        SINO3= SIN(UT3)
        COSO3= COS(UT3)
        TANO3= SINO3/ COSO3
        XB3 = RM3(KMAG)* COS(UMEG3) - RM
        YB3 = RM3(KMAG)* SIN(UMEG3)
        XA3 = U23(KMAG) * COSO3 + XB3
        YA3 = YB3 + U23(KMAG) * SINO3
        XCC = U13(KMAG) * COSO3 + XB3
        YCC = YB3 + U13(KMAG) * SINO3
        XD3 = XCC + R13(KMAG) * SINO3
        YD3 = YCC - R13(KMAG) * COSO3
        XE3 = XA3 + R23(KMAG) * SINO3
        YE3 = YA3 - R23(KMAG) * COSO3
        SIN23 = SINO3**2
        COS23 = COSO3**2
        SICO3 = SINO3 * COSO3
      ENDIF
C---------- endif third EFB

C----- Compute field at X, Y
C      Coordinates of current position

      DO 1 JRO = 1,NN
        ROJ = RO + DRO * DBLE(NN-JRO-INT(NN/2))

        DO 1 ITTA = 1,NN
          TTAI = TTA - DTTA * DBLE(NN-ITTA-INT(NN/2))
          ZETA = ACN(KMAG) - TTAI

C          X = ROJ * COS(ZETA) - RM
          X = ROJ * COS(ZETA) - RRM(KMAG)
          Y = ROJ * SIN(ZETA)

C         ... AX (CX) = COTE X DE LA PROJECTION DE A (C) SUR LA DROITE (MX)
C         (M=POINT COURANT), ORTHOGONALEMENT A LA DROITE (ABC)
          AX = ((YA - Y) * TANO ) + XA
          CX = ((YC - Y) * TANO ) + XC
          AXS = ((YAS - Y) * TANOS ) + XAS
          CXS = ((YCS - Y) * TANOS ) + XCS
          IF(NBFACE(KMAG) .EQ. 3) THEN
            AX3 = ((YA3 - Y) * TANO3 ) + XA3
            CX3 = ((YCC - Y) * TANO3 ) + XCC
          ENDIF
C
C POSITION DU POINT (X,Y) / FACE ENTREE
          IF     ( X.LT.CX .AND. X.LT.AX )  THEN
C           REGION DE COURBURE R1
            RR = R1(KMAG)
            XO = XD
            YO = YD
            R = ABS(RR)
            D = R - SQRT((X - XO)**2 + (Y - YO)**2 )
            D =( -D * RR/R + SHIFTE(KMAG) )
          ELSE IF( X.GE.CX .AND. X.LE.AX )  THEN
C            REGION LINEAIRE
            XO = SICO * (Y - YB) + XB * SIN2 + X * COS2
            YO = SICO * (X - XB) + YB * COS2 + Y * SIN2
            YL = YB + ((X - XB) * TANO )
C           ... YL = COTE Y DE LA PROJECTION DE B SUR LA DROITE (MY)
C           (M=POINT COURANT), PARALLELEMENT A LA DROITE (ABC)
            D = SQRT((X - XO)**2 + (Y - YO)**2 )
            IF( Y .LE. YL .OR. D .LE. 1.D-6 ) D = -D
            D=( D + SHIFTE(KMAG) )
          ELSE IF( X.GT.CX .AND. X.GT.AX )  THEN
C           REGION DE COURBURE R2
            RR = R2(KMAG)
            XO = XE
            YO = YE
            R = ABS(RR)
            D = R - SQRT((X - XO)**2 + (Y - YO)**2 )
            D=( -D * RR/R + SHIFTE(KMAG) )
          ELSE
C           ERREUR  DE  DONNEES  FACE  ENTREE
            IF(NRES.GT.0) WRITE(NRES,104)
 104        FORMAT(/,5X,10('*'),
     >      'Pgm dipsi. ERREUR PARAMETRES FACE ENTREE',/)
            IF(NRES.GT.0) WRITE(NRES,134) X,AX,CX,XA,YA,XC,YC
            GOTO  99
          ENDIF

          GAP =  LAMBDE(KMAG) * (RM/RO)**QSIE(KMAG)
          IF(GAP .EQ. 0.D0) THEN
            IF(D.LE.0.D0) THEN
              FE=1.D0
            ELSE
              FE=0.D0
            ENDIF
          ELSE
            D = D/GAP
            P=CE(KMAG,1)+(CE(KMAG,2)+(CE(KMAG,3)+(CE(KMAG,4)+
     >                       (CE(KMAG,5)+CE(KMAG,6)*D)*D)*D)*D)*D
            IF    (P .GE.  PLIM) THEN
              FE = 0.D0
            ELSEIF(P .LE. -PLIM) THEN
              FE = 1.D0
            ELSE
              FE = 1.D0/(1.D0+EXP(P))
            ENDIF
          ENDIF

C POSITION DU POINT (X,Y) / FACE SORTIE
          IF     ( X.LT.CXS .AND. X.LT.AXS )  THEN
C           ... REGION DE COURBURE R1
            RRS = R1S(KMAG)
            XO = XDS
            YO = YDS
            R = ABS(RRS)
            D = R - SQRT((X - XO)**2 + (Y - YO)**2 )
            D=( D * RRS/R + SHIFTS(KMAG) )
          ELSE IF( X.GE.CXS .AND. X.LE.AXS )  THEN
C           ... REGION LINEAIRE
            XO = SICOS * (Y - YBS) + XBS * SIN2S + X * COS2S
            YO = SICOS * (X - XBS) + YBS * COS2S + Y * SIN2S
            YL = YBS + ((X - XBS) * TANOS )
            D = SQRT((X - XO)**2 + (Y - YO)**2 )
            IF( Y .GE. YL .OR. D .LE.1.D-6 ) D = -D
            D=( D + SHIFTS(KMAG) )
          ELSE IF( X.GT.CXS .AND. X.GT.AXS )  THEN
C           ... REGION DE COURBURE R2
            RRS = R2S(KMAG)
            XO = XES
            YO = YES
            R = ABS(RRS)
            D = R - SQRT((X - XO)**2 + (Y - YO)**2 )
            D=( D * RRS/R + SHIFTS(KMAG) )
          ELSE
            IF(NRES.GT.0) WRITE(NRES,114)
114         FORMAT(/,5X,10('*'),
     >      'Pgm dipsi. ERREUR PARAMETRES FACE SORTIE',/)
            IF(NRES.GT.0) WRITE(NRES,134) X,AXS,CXS,XAS,YAS,XCS,YCS
            GOTO  99
          ENDIF

          GAP =  LAMBDS(KMAG) * (RM/RO)**QSIS(KMAG)
          IF(GAP .EQ. 0.D0) THEN
            IF(D.LE.0.D0) THEN
              FS=1.D0
            ELSE
              FS=0.D0
            ENDIF
          ELSE
            D = D/GAP
            P=CS(KMAG,1)+(CS(KMAG,2)+(CS(KMAG,3)+(CS(KMAG,4)+
     >                        (CS(KMAG,5)+CS(KMAG,6)*D)*D)*D)*D)*D
            IF    (P .GE.  PLIM) THEN
              FS = 0.D0
            ELSEIF(P .LE. -PLIM) THEN
              FS = 1.D0
            ELSE
              FS = 1.D0/(1.D0+EXP(P))
            ENDIF
          ENDIF

C POSITION DU POINT (X,Y) / FACE EXTERNE (SEULEMENT SI OPTION 3 FACES)
          IF(NBFACE(KMAG) .EQ. 3) THEN
            IF     ( X.LT.CX3 .AND. X.LT.AX3 )  THEN
C             ... REGION DE COURBURE R1
              RR3 = R13(KMAG)
              XO = XD3
              YO = YD3
              R = ABS(RR3)
              D = R - SQRT((X - XO)**2 + (Y - YO)**2 )
              D=( D * RR3/R + SHIFT3(KMAG) )
            ELSE IF( X.GE.CX3 .AND. X.LT.AX3 )  THEN
C             ... REGION LINEAIRE
              XO = SICO3 * (Y - YB3) + XB3 * SIN23 + X * COS23
              YO = SICO3 * (X - XB3) + YB3 * COS23 + Y * SIN23
              YL = YB3 + ((X - XB3) * TANO3 )
              D = SQRT((X - XO)**2 + (Y - YO)**2 )
              IF( Y .GE. YL .OR. D .LE. 1.D-6 ) D = -D
              D=( D + SHIFT3(KMAG) )
            ELSE IF( X.GE.CX3 .AND. X.GE.AX3 )  THEN
C             ... REGION DE COURBURE R2
              RR3 = R23(KMAG)
              XO = XE3
              YO = YE3
              R = ABS(RR3)
              D = R - SQRT((X - XO)**2 + (Y - YO)**2 )
              D=( D * RR3/R + SHIFT3(KMAG) )
            ELSE
              IF(NRES.GT.0) THEN
                WRITE(NRES,124)
124             FORMAT(/,5X,9('*'),
     >          'Pgm dipsi. ERREUR PARAMETRES FACE LATERALE',/)
                WRITE(NRES,134) X,AX3,CX3,XA3,YA3,XCC,YCC
134             FORMAT(/,5X,7G12.4)
              ENDIF
              GOTO  99
            ENDIF

            GAP =  LAMBD3(KMAG)
            IF(GAP .EQ. 0.D0) THEN
              IF(D.LE.0.D0) THEN
                F3=1.D0
              ELSE
                F3=0.D0
              ENDIF
            ELSE
              D = D/GAP
              P=CC(KMAG,1)+(CC(KMAG,2)+(CC(KMAG,3)+(CC(KMAG,4)+
     >                           (CC(KMAG,5)+CC(KMAG,6)*D)*D)*D)*D)*D
              IF    (P .GE.  PLIM) THEN
                F3 = 0.D0
              ELSEIF(P .LE. -PLIM) THEN
                F3 = 1.D0
              ELSE
                F3 = 1.D0/(1.D0+EXP(P))
              ENDIF
            ENDIF
          ENDIF

          F = FE * FS
          IF(NBFACE(KMAG) .EQ. 3) F = F * F3

          ROI= ROJ-RRM(KMAG)

          IF    (ITYPF.EQ.0) THEN
C Default
            ROI= ROI/RRM(KMAG)
            SF = 1.D0
            IF(IND(KMAG).NE.0) THEN
              PROI = 1.D0
              DO KND = 1, IND(KMAG)
                PROI = PROI * ROI
                SF = SF + CIND(KMAG,KND) * PROI
              ENDDO
            ENDIF
            FTAB(ITTA,JRO) = FTAB(ITTA,JRO) + F*HNORM(KMAG) * SF
          ELSEIF(ITYPF.EQ.1) THEN
C B0 can be zero
            SF = HNORM(KMAG)
            IF(IND(KMAG).NE.0) THEN
              PROI = 1.D-1 ! convert CIND from T/m^* to kG/cm^*
              DO KND = 1, IND(KMAG)
                PROI = PROI * ROI
                SF = SF + CIND(KMAG,KND) * PROI
              ENDDO
            ENDIF
            FTAB(ITTA,JRO) = FTAB(ITTA,JRO) + F * SF
          ELSEIF(ITYPF.GE.7) THEN
C Sandro's type of field for AGS upgrade FFAG
            ROI = ROI * 1.D-2
            ROI2 = ROI  * ROI
            ROI3 = ROI2 * ROI
            ROI4 = ROI3 * ROI
            ROI5 = ROI4 * ROI
            SF = HNORM(KMAG) + CIND(KMAG,1) * ROI

            IF    (ITYPF.EQ.7) THEN
C D-type dipole. TTRF is at half the sector angle.

              TTRF = ACN(KMAG)-0.5D0*(UMEG+UMEGS)
              TTA2 = (TTAI-TTRF) * (TTAI-TTRF) * 1.D6
              TTA4 = TTA2 * TTA2
              TTA6 = TTA4 * TTA2
              TTA8 = TTA6 * TTA2
              sf =  sf   +
     >        ROI2*(ad20 +ad22*TTA2 +ad24*TTA4 +ad26*TTA6 +ad28*TTA8) +
     >        ROI3*(ad30 +ad32*TTA2 +ad34*TTA4 +ad36*TTA6 +ad38*TTA8) +
     >        ROI4*(ad40 +ad42*TTA2 +ad44*TTA4 +ad46*TTA6 +ad48*TTA8) +
     >        ROI5*(ad50 +ad52*TTA2 +ad54*TTA4 +ad56*TTA6 +ad58*TTA8)
              IF(IND(KMAG).LT.0)  SF = -SF

            ELSEIF(ITYPF.EQ.8 .OR. ITYPF.EQ.9) THEN
C F-type dipole
              IF    (ITYPF.EQ.8) THEN
C First type F-type dipole
                TTRF = ACN(KMAG)-UMEGS
              ELSEIF(ITYPF.EQ.9) THEN
C Second type F-type dipole
                TTRF = ACN(KMAG)-UMEG
              ENDIF

              TTA2 = (TTAI-TTRF)**2  * 1.D6
              TTA4 = TTA2 * TTA2
              TTA6 = TTA4 * TTA2
              TTA8 = TTA6 * TTA2
              SF =  SF   +
     >        ROI2*(AF20 +AF22*TTA2 +AF24*TTA4 +AF26*TTA6 +AF28*TTA8) +
     >        ROI3*(AF30 +AF32*TTA2 +AF34*TTA4 +AF36*TTA6 +AF38*TTA8) +
     >        ROI4*(AF40 +AF42*TTA2 +AF44*TTA4 +AF46*TTA6 +AF48*TTA8) +
     >        ROI5*(AF50 +AF52*TTA2 +AF54*TTA4 +AF56*TTA6 +AF58*TTA8)
              IF(IND(KMAG).LT.0)  SF = -SF

            ENDIF

            FTAB(ITTA,JRO) = FTAB(ITTA,JRO) + F * SF
          ENDIF

 1    CONTINUE
C--------- end loop on   NN x tta  &  NN x ro

      IF(KMAG.LT.NBMAG) GOTO 20
C------------------------------
      RETURN

C-------------------------------------------------------------------
C----- Analytic calculation of B & derivatives ---------------------
      ENTRY DIPSFA(IDB,TTA,RO,
     >                        BZ0)
      KMAG = 0
 30   CONTINUE
      KMAG = KMAG+1

      QGAP=QSIE(KMAG)
      QGAPS=QSIS(KMAG)
      HO=HNORM(KMAG)
      IF    (ITYPF.EQ.0) THEN
        COEF1 = CIND(KMAG,1)
        COEF2 = CIND(KMAG,2)
        COEF3 = CIND(KMAG,3)
        COEF4 = CIND(KMAG,4)
        COEF5 = CIND(KMAG,5)
      ELSEIF(ITYPF.EQ.1) THEN
        COEF1 = CIND(KMAG,1)*RM    *1.D-1
        COEF2 = CIND(KMAG,2)*RM**2 *1.D-1
        COEF3 = CIND(KMAG,3)*RM**3 *1.D-1
        COEF4 = CIND(KMAG,4)*RM**4 *1.D-1
        COEF5 = CIND(KMAG,5)*RM**5 *1.D-1
      ENDIF

      RR1 = R1(KMAG)
      RR2 = R2(KMAG)
      RR1S = R1S(KMAG)
      RR2S = R2S(KMAG)

C************************ FACE D'ENTREE *************************

C     ** POUR LES BESOINS DE LA SIMPLE PRECISION :
      IF(RR1*RR1  .GE. 1.D6) U1(KMAG) = -1.D6
      IF(RR2*RR2  .GE. 1.D6) U2(KMAG) =  1.D6

      UMEG = UMEGA(KMAG) * RAD
      TETA = THETA(KMAG) * RAD
      UT = UMEG - TETA
      SINO = SIN(UT)
      COSO = COS(UT)
      TANO = SINO / COSO

C     PROJECTIONS DU RAYON MOYEN SUR LES AXES X,Y
      XB = RM * ( COS(UMEG) - 1.D0)
      YB = RM * SIN(UMEG)
C     COORDONNEES DE L'EXTREMITE A DE LA PARTIE LINEAIRE DE LONG. U2(KMAG)
      XA = U2(KMAG) * COSO + XB
      YA = YB + U2(KMAG) * SINO
C     COORDONNEES DE L'EXTREMITE C DE LA PARTIE LINEAIRE DE LONG. U1
      XC = U1(KMAG) * COSO + XB
      YC = YB + U1(KMAG) * SINO
C     COORDONNEES DU CENTRE DE COURBURE DE RAYON RR1
      XD = XC + RR1 * SINO
      YD = YC - RR1 * COSO
C     COORDONNEES DU CENTRE DE COURBURE DE RAYON RR2
      XE = XA + RR2 * SINO
      YE = YA - RR2 * COSO

      SIN2 = SINO**2
      COS2 = COSO**2
      SICO = SINO * COSO

C************************ FACE DE SORTIE *********************
C     ** POUR LES BESOINS DE LA SIMPLE PRECISION :
      IF(RR1S*RR1S .GE. 1.D10) U1S(KMAG) = -1.D6
      IF(RR2S*RR2S .GE. 1.D10) U2S(KMAG) =  1.D6

      UMEGS = UMEGAS(KMAG) * RAD
      TETAS = THETAS(KMAG) * RAD
      UTS= UMEGS- TETAS
      SINOS= SIN(UTS)
      COSOS= COS(UTS)
      TANOS= SINOS/ COSOS

      XBS = RM * ( COS(UMEGS) - 1.D0)
      YBS = RM * SIN(UMEGS)
      XAS = U2S(KMAG) * COSOS + XBS
      YAS = YBS + U2S(KMAG) * SINOS
      XCS = U1S(KMAG) * COSOS + XBS
      YCS = YBS + U1S(KMAG) * SINOS
      XDS = XCS + RR1S * SINOS
      YDS = YCS - RR1S * COSOS
      XES = XAS + RR2S * SINOS
      YES = YAS - RR2S * COSOS

      SIN2S = SINOS**2
      COS2S = COSOS**2
      SICOS = SINOS * COSOS

C      Projections du point M sur les axes X et Y
       ZETA = ACN(KMAG)-TTA
       COSZ = COS(ZETA)
       SINZ = SIN(ZETA)
       X = RO*COSZ - RM
       Y = RO*SINZ

C         ... AX (CX) = COTE X DE LA PROJECTION DE A (C) SUR LA DROITE (MX)
C         (M=POINT COURANT), ORTHOGONALEMENT A LA DROITE (ABC)

      AX = ((YA - Y) * TANO ) + XA
      CX = ((YC - Y) * TANO ) + XC
      AXS = ((YAS - Y) * TANOS ) + XAS
      CXS = ((YCS - Y) * TANOS ) + XCS

C*****************   POSITION DE  M(X,Y)/ FACE D'ENTREE  *******************

C     On se place dans la zone lineaire i.e IF( X.GE.CX .AND. X.LE.AX )
C     Dans ce cas XO et YO ne sont pas des constantes

       XO = SICO * (Y - YB) + XB * SIN2 + X * COS2
       YO = SICO * (X - XB) + YB * COS2 + Y * SIN2
       YL = YB + ((X - XB) * TANO )
C       ... YL = COTE Y DE LA PROJECTION DE B SUR LA DROITE (MY)
C           (M=POINT COURANT), PARALLELEMENT A LA DROITE (ABC)
       D = SQRT((X - XO)**2 + (Y - YO)**2 )
       D2 = D * D
       ROOT   = D
       IF(ROOT.EQ.0.D0)
     > CALL ENDJOB('ERROR : root=0 in dipsfa. Try D=1e-99 in pgm',-99)
       SIGN = 1.D0
       IF( Y .LE. YL .OR. D .LE. 1.D-6 ) THEN
          D=-D
          SIGN =- 1.D0
       ENDIF
C       D=( D + SHIFTE )

C**************** CALCUL DES DERIVEES POUR LA FACE D'ENTREE *************

       GAPO = LAMBDE(KMAG)

       C0 = CE(KMAG,1)
       C1 = CE(KMAG,2)
       C2 = CE(KMAG,3)
       C3 = CE(KMAG,4)
       C4 = CE(KMAG,5)
       C5 = CE(KMAG,6)

       XRO    =  COSZ
       X2RO   =  0.D0
       XTTA   =  RO*SINZ
       X2TTA  = -RO*COSZ
       XRTA   =  SINZ

       YRO    =  SINZ
       Y2RO   =  0.D0
       YTTA   = -RO*COSZ
       Y2TTA  = -RO*SINZ
       YRTA   = -COSZ

       XORO   = SICO*YRO     + COS2*XRO
       XO2RO  = SICO*Y2RO    + COS2*X2RO
       XOTTA  = SICO*YTTA    + COS2*XTTA
       XO2TA  = SICO*Y2TTA   + COS2*X2TTA
       XORTA  = SICO*YRTA    + COS2*XRTA

       YORO   = SICO*XRO     + SIN2*YRO
       YO2RO  = SICO*X2RO    + SIN2*Y2RO
       YOTTA  = SICO*XTTA    + SIN2*YTTA
       YO2TA  = SICO*X2TTA   + SIN2*Y2TTA
       YORTA  = SICO*XRTA    + SIN2*YRTA

       DX     = X      -  XO
       DXRO   = XRO    -  XORO
       DX2RO  = X2RO   -  XO2RO
       DXTTA  = XTTA   -  XOTTA
       DX2TA  = X2TTA  -  XO2TA
       DXRTA  = XRTA   -  XORTA

       DY     = Y      -  YO
       DYRO   = YRO    -  YORO
       DY2RO  = Y2RO   -  YO2RO
       DYTTA  = YTTA   -  YOTTA
       DY2TA  = Y2TTA  -  YO2TA
       DYRTA  = YRTA   -  YORTA

       ROOT2  = ROOT*ROOT
       ROOT3  = ROOT2*ROOT

       D1RO   =   (DXRO*DX  + DYRO*DY)  / ROOT
       D2RO   = - (DX*DXRO  + DY*DYRO)**2 / ROOT3
     >      + (DXRO*DXRO + DYRO*DYRO + DX*DX2RO + DY*DY2RO) / ROOT
       D1TTA  =   (DXTTA*DX + DYTTA*DY)  /ROOT
       D2TTA  = - (DX*DXTTA + DY*DYTTA)**2 /ROOT3
     >      + (DXTTA*DXTTA + DYTTA*DYTTA + DX*DX2TA + DY*DY2TA)
     >      / ROOT
       DRTA  = - (DX*DXTTA + DY*DYTTA)*(DX*DXRO + DY*DYRO)
     >      /  ROOT3
     >      + (DXTTA*DXRO + DYTTA*DYRO + DX*DXRTA + DY*DYRTA)
     >      /  ROOT

       GAP =  GAPO * (RM/RO)**QGAP
       GAP2 = GAP*GAP
       GAP3 = GAP2*GAP

       DL=D/GAP
       UL=1.D0/GAP
       UL2 = UL*UL
       P     = C0+(C1+(C2+(C3+(C4+C5*DL)*DL)*DL)*DL)*DL
       PP    = C1+(2*C2+(3*C3+(4*C4+5*C5*DL)*DL)*DL)*DL
       PPP   = 2.D0*C2+(6.D0*C3+(12.D0*C4+20.D0*C5*DL)*DL)*DL
       EP=EXP(P)

       GARO  = -GAPO*(QGAP/RM)*(RM/RO)**(QGAP+1)
       GA2RO =  GAPO*(QGAP*(QGAP+1.D0)/(RM*RM))*(RM/RO)**(QGAP+2.D0)
       DD  =  ((SIGN*D1RO)*GAP-D*GARO)/ GAP2
       D2D = ((SIGN*D2RO)*GAP2-D*GAP*GA2RO
     > +2.D0*GARO*GARO*D-2.D0*GAP*(SIGN*D1RO)*GARO)/GAP3
       PPTA  = UL*(SIGN*D1TTA)*PPP
       PPR   = DD*PPP

       PRO    =  DD * PP
       P2RO   =  D2D*PP + DD*DD*PPP
       PTTA  = SIGN * UL*D1TTA * PP
       P2TTA = UL*((SIGN*D2TTA)*PP + (SIGN*D1TTA)*PPTA)
       PRTA   = (-UL2*GARO*(SIGN*D1TTA) + UL*(SIGN*DRTA))*PP
     >         +UL*(SIGN*D1TTA)*PPR

       IF(GAP .EQ. 0.D0) THEN
         FRO    = 0.D0
         F2RO   = 0.D0
         FTTA   = 0.D0
         F2TTA  = 0.D0
         FRTA   = 0.D0
         IF(D.LE.0.D0) THEN
           F=1.D0
         ELSE
           F=0.D0
         ENDIF
       ELSE
         IF  (P .GE. PLIM) THEN
           F      = 0.D0
           FRO    = 0.D0
           F2RO   = 0.D0
           FTTA   = 0.D0
           F2TTA  = 0.D0
           FRTA   = 0.D0
         ELSEIF(P .LE. -PLIM) THEN
           F      = 1.D0
           FRO    = 0.D0
           F2RO   = 0.D0
           FTTA   = 0.D0
           F2TTA  = 0.D0
           FRTA   = 0.D0
         ELSE
           UEP  = 1.D0+EP
           UMEP = 1.D0-EP
           F  = 1.D0/UEP
           F2 = F  * F
           F3 = F2 * F
           PRO2  = PRO*PRO
           PTTA2 = PTTA*PTTA
           FRO   = -PRO*EP * F2
           F2RO  = -(EP * F3)*(P2RO*UEP+PRO2*UMEP)
           FTTA  = -PTTA*EP * F2
           F2TTA = -(EP * F3)*(P2TTA*UEP+PTTA2*UMEP)
           FRTA  = -(EP * F3)*(PRTA*UEP+(PTTA*PRO*UMEP))
         ENDIF
       ENDIF

C       IF (KIRD.EQ.1) THEN
       IF (IDB.GE.4) THEN
          X3RO   =  0.D0
          X4RO   =  0.D0
          X3TTA  = -RO*SINZ
          X4TTA  =  RO*COSZ
          X2RTA  =  0.D0
          XR2TA  = -COSZ
          X3RTA  =  0.D0
          XR3TA  = -SINZ
          X2R2TA =  0.D0

          Y3RO   =  0.D0
          Y4RO   =  0.D0
          Y3TTA  =  RO*COSZ
          Y4TTA  =  RO*SINZ
          Y2RTA  =  0.D0
          YR2TA  = -SINZ
          Y3RTA  =  0.D0
          YR3TA  =  COSZ
          Y2R2TA =  0.D0

          XO3RO   = SICO*Y3RO    + COS2*X3RO
          XO4RO   = SICO*Y4RO    + COS2*X4RO
          XO3TA   = SICO*Y3TTA   + COS2*X3TTA
          XO4TA   = SICO*Y4TTA   + COS2*X4TTA
          XO2RTA  = SICO*Y2RTA   + COS2*X2RTA
          XOR2TA  = SICO*YR2TA   + COS2*XR2TA
          XO3RTA  = SICO*Y3RTA   + COS2*X3RTA
          XOR3TA  = SICO*YR3TA   + COS2*XR3TA
          XO2R2TA = SICO*Y2R2TA  + COS2*X2R2TA

          YO3RO   = SICO*X3RO    + SIN2*Y3RO
          YO4RO   = SICO*X4RO    + SIN2*Y4RO
          YO3TA   = SICO*X3TTA   + SIN2*Y3TTA
          YO4TA   = SICO*X4TTA   + SIN2*Y4TTA
          YO2RTA  = SICO*X2RTA   + SIN2*Y2RTA
          YOR2TA  = SICO*XR2TA   + SIN2*YR2TA
          YO3RTA  = SICO*X3RTA   + SIN2*Y3RTA
          YOR3TA  = SICO*XR3TA   + SIN2*YR3TA
          YO2R2TA = SICO*X2R2TA  + SIN2*Y2R2TA

          DX3RO = X3RO   -  XO3RO
          DX4RO = X4RO   -  XO4RO
          DX3TA = X3TTA  -  XO3TA
          DX4TA = X4TTA  -  XO4TA

          DY3RO = Y3RO   -  YO3RO
          DY4RO = Y4RO   -  YO4RO
          DY3TA = Y3TTA  -  YO3TA
          DY4TA = Y4TTA  -  YO4TA

          DX2RTA  =  X2RTA  -  XO2RTA
          DXR2TA  =  XR2TA  -  XOR2TA
          DX3RTA  =  X3RTA  -  XO3RTA
          DXR3TA  =  XR3TA  -  XOR3TA
          DX2R2TA =  X2R2TA -  XO2R2TA

          DY2RTA  =  Y2RTA  -  YO2RTA
          DYR2TA  =  YR2TA  -  YOR2TA
          DY3RTA  =  Y3RTA  -  YO3RTA
          DYR3TA  =  YR3TA  -  YOR3TA
          DY2R2TA =  Y2R2TA -  YO2R2TA

          ROOT5  = ROOT2*ROOT3
          ROOT7  = ROOT5*ROOT2

          D3RO   =   (3.D0*(DX*DXRO + DY*DYRO)**3) / ROOT5
     >     - (3.D0*(DX*DXRO + DY*DYRO)*(DXRO**2 + DYRO**2 + DX*DX2RO
     >     +  DY*DY2RO)) / ROOT3
     >     + (3.D0*DXRO*DX2RO + 3.D0*DYRO*DY2RO + DX*DX3RO + DY*DY3RO)
     >     /  ROOT
          D4RO   =   -15.D0*(DX* DXRO + DY* DYRO)**4 / ROOT7
     >     + 9.D0*(DX*DXRO + DY*DYRO)**2.D0*(2.D0*DXRO**2 + 2.D0*DYRO**2
     >     + 2.D0*DX*DX2RO + 2.D0*DY*DY2RO) / ROOT5
     >     - 3.D0*(DXRO**2 + DYRO**2 + DX*DX2RO + DY*DY2RO)**2
     >     / ROOT3
     >     - 4.D0*(DX*DXRO + DY*DYRO)*(3.D0*DXRO*DX2RO + 3.D0*DYRO*DY2RO
     >     + DX*DX3RO + DY*DY3RO) / ROOT3
     >     + (3.D0*DX2RO**2 + 3.D0*DY2RO**2 + 4.D0*DXRO*DX3RO
     >     + 4.D0*DYRO*DY3RO + DX*DX4RO + DY*DY4RO) / ROOT
          D3TTA  =   (3.D0*(DX*DXTTA + DY*DYTTA)**3) / ROOT5
     >     - (3.D0*(DX*DXTTA + DY*DYTTA)*(DXTTA**2 + DYTTA**2
     >     + DX*DX2TA + DY*DY2TA)) / ROOT3
     >     + (3.D0*DXTTA*DX2TA + 3.D0*DYTTA*DY2TA+DX*DX3TA+DY*DY3TA)
     >     /  ROOT
          D4TTA   =   -15*(DX* DXTTA + DY* DYTTA)**4 / ROOT7
     >     + 9.D0*(DX*DXTTA+ DY*DYTTA)**2.D0*(2*DXTTA**2+ 2.D0*DYTTA**2
     >     + 2.D0*DX*DX2TA + 2.D0*DY*DY2TA) / ROOT5
     >     - 3.D0*(DXTTA**2 + DYTTA**2 + DX*DX2TA + DY*DY2TA)**2
     >     / ROOT3
     >     - 4.D0*(DX*DXTTA+DY*DYTTA)*(3.D0*DXTTA*DX2TA+3.D0*DYTTA*DY2TA
     >     + DX*DX3TA + DY*DY3TA) / ROOT3
     >     + (3.D0*DX2TA**2 + 3.D0*DY2TA**2 + 4.D0*DXTTA*DX3TA
     >     + 4.D0*DYTTA*DY3TA + DX*DX4TA + DY*DY4TA) / ROOT
          D2RTA =   3.D0*(DX*DXTTA + DY*DYTTA)*(DX*DXRO + DY*DYRO)**2
     >     /  ROOT5
     >     -  (DX*DXRO + DY*DYRO)*(2.D0*DXTTA*DXRO + 2.D0*DYTTA*DYRO
     >     +  2.D0*DX*DXRTA + 2.D0*DY*DYRTA) / ROOT3
     >     -  ((DX*DXTTA + DY*DYTTA)*(DXRO**2 + DYRO**2 + DX*DX2RO
     >     +  DY*DY2RO)) / ROOT3
     >     +  (2.D0*DXRO*DXRTA + 2.D0*DYRO*DYRTA + DXTTA*DX2RO
     >     +  DYTTA*DY2RO +  DX*DX2RTA + DY*DY2RTA) / ROOT
          DR2TA = -(DX*DXTTA + DY*DYTTA)*(2.D0*DX*DXRTA+ 2.D0*DXRO*DXTTA
     >     +  2.D0*DY*DYRTA + 2.D0*DYRO*DYTTA) / ROOT3
     >     +  (DX*DXRO + DY*DYRO)*(3.D0*(DX*DXTTA + DY*DYTTA)**2
     >     / ROOT5
     >         - (DX*DX2TA+DXTTA**2+DY*DY2TA+DYTTA**2) / ROOT3)
     >         + (DX*DXR2TA + DX2TA*DXRO + 2.D0*DXRTA*DXTTA + DY*DYR2TA
     >         + DY2TA*DYRO + 2.D0*DYRTA*DYTTA) / ROOT
          D3RTA = - 15.D0*(DX*DXRO + DY*DYRO)**3 * (DX*DXTTA + DY*DYTTA)
     >         / ROOT7
     >         - (DX*DX3RO + 3.D0*DX2RO*DXRO+ DY*DY3RO+ 3.D0*DY2RO*DYRO)
     >         * (DX*DXTTA + DY*DYTTA) / ROOT3
     >         + 9.D0*(DX*DXRO + DY*DYRO)*(DX*DX2RO + DXRO**2 + DY*DY2RO
     >         + DYRO**2)*(DX*DXTTA + DY*DYTTA)/ ROOT5
     >         - 3.D0*(DX*DXRO + DY*DYRO)*(DX*DX2RTA + 2.D0*DXRO*DXRTA
     >         +  DX2RO*DXTTA + DY*DY2RTA + 2.D0*DYRO*DYRTA+DY2RO*DYTTA)
     >         /  ROOT3
     >         + (DX*DX3RTA + 3*DX2RTA*DXRO+3*DX2RO*DXRTA+DX3RO*DXTTA
     >         +  DY*DY3RTA + 3*DY2RTA*DYRO+3*DY2RO*DYRTA+DY3RO*DYTTA)
     >         /  ROOT
     >         +  9.D0*(DX*DXRO + DY*DYRO)**2 * (DX*DXRTA + DXRO*DXTTA
     >         +  DY*DYRTA + DYRO*DYTTA) / ROOT5
     >         -  3.D0*(DX*DX2RO + DXRO**2 +DY*DY2RO+ DYRO**2)*(DX*DXRTA
     >         +  DXRO*DXTTA + DY*DYRTA + DYRO*DYTTA) / ROOT3
          DR3TA = (DX*DXR3TA + DX3TA*DXRO+3*DX2TA*DXRTA+3*DXR2TA*DXTTA
     >       + DY*DYR3TA + DY3TA*DYRO+3*DY2TA*DYRTA+3*DYR2TA*DYTTA)
     >       / ROOT
     >       - 3.D0*(DX*DXTTA + DY*DYTTA)*(DX*DXR2TA + DX2TA*DXRO
     >       + 2.D0*DXRTA*DXTTA + DY*DYR2TA + DY2TA*DYRO+2*DYRTA*DYTTA)
     >       / ROOT3
     >       + 3.D0*(DX*DXRTA + DXRO*DXTTA + DY*DYRTA + DYRO*DYTTA)
     >       * ((3.D0*(DX*DXTTA + DY*DYTTA)**2)/ ROOT5
     >       - (DX*DX2TA + DXTTA**2 + DY*DY2TA + DYTTA**2)
     >       / ROOT3)
     >       + (DX*DXRO + DY*DYRO)*((-15.D0*(DX*DXTTA + DY*DYTTA)**3)
     >       / ROOT7
     >       - (DX*DX3TA+ 3.D0*DX2TA*DXTTA+ DY*DY3TA+ 3.D0*DY2TA*DYTTA)
     >       /  ROOT3
     >       + (9.D0*(DX*DXTTA + DY*DYTTA)*(DX*DX2TA+DXTTA**2+DY*DY2TA
     >       + DYTTA**2))/ROOT5)
          D2R2TA =-(DX*DXTTA+ DY*DYTTA)*(2.D0*DX*DX2RTA+ 4.D0*DXRO*DXRTA
     >         + 2*DX2RO*DXTTA  + 2*DY*DY2RTA + 4*DYRO*DYRTA
     >         + 2*DY2RO*DYTTA) / ROOT3
          D2R2TA = D2R2TA
     >         + 3*(DX*DXRO + DY*DYRO)*(2*DX*DXTTA + 2*DY*DYTTA)
     >         * (2*DX*DXRTA + 2*DXRO*DXTTA + 2*DY*DYRTA
     >         + 2*DYRO*DYTTA) / ROOT5
          D2R2TA = D2R2TA
     >         - (2*(DX*DXRTA + DXRO*DXTTA + DY*DYRTA+DYRO*DYTTA)**2
     >         + 2*(DX*DXRO + DY*DYRO)*(DX*DXR2TA + DX2TA*DXRO
     >         + 2*DXRTA*DXTTA + DY*DYR2TA + DY2TA*DYRO
     >         + 2*DYRTA*DYTTA))/ROOT3
          D2R2TA = D2R2TA
     >         + (DX*DXRO + DY*DYRO)**2*(-15*(DX*DXTTA +DY*DYTTA)**2
     >         / ROOT7
     >         + 3*(DX*DX2TA + DXTTA**2 + DY*DY2TA
     >         + DYTTA**2)/ROOT5)
     >         + (DX*DX2RO + DXRO**2 + DY*DY2RO + DYRO**2)
     >         * (3*(DX*DXTTA + DY*DYTTA)**2/ROOT5
     >         - (DX*DX2TA + DXTTA**2 + DY*DY2TA + DYTTA**2)
     >         /  ROOT3)
     >         +  (DX*DX2R2TA +DX2RO*DX2TA+2*DXR2TA*DXRO+2*DXRTA**2
     >         +  2*DX2RTA*DXTTA + DY*DY2R2TA + DY2RO*DY2TA
     >         +  2*DYR2TA*DYRO + 2*DYRTA**2 + 2*DY2RTA*DYTTA)
     >         /  ROOT

          UL3 = UL2*UL
          UL4 = UL3*UL
          PPPP  = 6.D0*C3+(24.D0*C4+60.D0*C5*DL)*DL
          PPPPP = 24.D0*C4+120.D0*C5*DL
          GAP4 = GAP3*GAP
          GAP5 = GAP4*GAP
          QGAP2 = QGAP*QGAP
          QGAP3 = QGAP2*QGAP

          GA3RO = -GAPO*QGAP*((2.D0+3.D0*QGAP+QGAP2)/(RM*RM*RM))
     >           *(RM/RO)**(QGAP+3.D0)
          GA4RO = GAPO*QGAP*((6.D0+11.D0*QGAP+6.D0*QGAP2+QGAP3)
     >         /(RM*RM*RM*RM))*(RM/RO)**(QGAP+4.D0)
          D3D= (6*GARO**2 *(GAP*(SIGN*D1RO)-D*GARO)
     >         -3.D0*GAP2*(GARO*(SIGN*D2RO)+(SIGN*D1RO)*GA2RO)
     >         +6.D0*D*GAP*GARO*GA2RO+GAP2*(GAP*(SIGN*D3RO)-D*GA3RO))
     >           /GAP4
          D4D= (6*GA2RO*GAP2*(D*GA2RO-(SIGN*D2RO)*GAP)
     >         -GAP3*(D*GA4RO-GAP*(SIGN*D4RO))
     >         +24*(SIGN*D1RO)*GA2RO*GAP2*GARO
     >         + 8*D*GA3RO*GAP2*GARO
     >         -12*GAP*GARO*GARO*(3*D*GA2RO-(SIGN*D2RO)*GAP)
     >         -4.D0*GAP3*((SIGN*D1RO)*GA3RO+(SIGN*D3RO)*GARO)
     >         -24.D0*GARO*GARO*GARO*((SIGN*D1RO)*GAP-D*GARO)) / GAP5
          PPPTA = UL*(SIGN*D1TTA)*PPPP
          PP2TA = UL*((SIGN*D2TTA)*PPP
     >         + (SIGN*D1TTA)*UL*(SIGN*D1TTA)*PPPP)
          PP3TA = UL*((SIGN*D3TTA)*PPP
     >         +3.D0*(SIGN*D2TTA)*UL*(SIGN*D1TTA)*PPPP
     >         + UL2*(SIGN*D1TTA)*(SIGN*D1TTA)*(SIGN*D1TTA)*PPPPP )
          PPPR  = DD*PPPP
          PP2R  = D2D*PPP + DD*DD*PPPP
          PPP2R = D2D*PPPP+ DD*DD*PPPPP
          PP3R  = D3D*PPP + 3*D2D*DD*PPPP + DD*DD*DD*PPPPP
          PPRTA = -UL2*GARO*(SIGN*D1TTA)*PPP
     >         +UL*(SIGN*DRTA)*PPP
     >         +UL*(SIGN*D1TTA)*DD*PPPP
          PPPRTA = -UL2*GARO*(SIGN*D1TTA)*PPPP
     >         +UL*(SIGN*DRTA)*PPPP
     >         +UL*(SIGN*D1TTA)*DD*PPPPP
          PP2RTA = (-UL2*2*(SIGN*DRTA)*GARO+UL*(SIGN*D2RTA)
     >         +(SIGN*D1TTA)*(2*UL3*GARO*GARO-UL2*GA2RO))*PPP
     >         + 2*(UL*(SIGN*DRTA)-UL2*(SIGN*D1TTA)*GARO)*PPPR
     >         + UL*(SIGN*D1TTA)*PPP2R
          PPR2TA  = (UL*(SIGN*DR2TA)-UL2*GARO*(SIGN*D2TTA))*PPP
     >         +(UL*(SIGN*DRTA)-UL2*(SIGN*D1TTA)*GARO)*PPPTA
     >         +(UL*(SIGN*D2TTA))*PPPR
     >         +UL*(SIGN*D1TTA)*PPPRTA

          P3RO   =  D3D*PP + 3*D2D*DD*PPP + DD*DD*DD*PPPP
          P4RO   =  D4D*PP + 4*D3D*DD*PPP + 3*D2D*D2D*PPP
     >         + 6*D2D*DD*DD*PPPP + DD*DD*DD*DD*PPPPP
          P3TTA = UL*((SIGN*D3TTA)*PP + 2*(SIGN*D2TTA)*PPTA
     >         + (SIGN*D1TTA)*PP2TA)
          P4TTA = UL*((SIGN*D4TTA)*PP + 3*(SIGN*D3TTA)*PPTA
     >         + 3*(SIGN*D2TTA)*PP2TA + (SIGN*D1TTA)*PP3TA)
          P2RTA  = (-UL2*2*(SIGN*DRTA)*GARO+UL*(SIGN*D2RTA)
     >         +(SIGN*D1TTA)*(2*UL3*GARO*GARO-UL2*GA2RO))*PP
     >         + 2*(UL*(SIGN*DRTA)-UL2*(SIGN*D1TTA)*GARO)*PPR
     >         + UL*(SIGN*D1TTA)*PP2R
          PR2TA  = (UL*(SIGN*DR2TA)-UL2*GARO*(SIGN*D2TTA))*PP
     >         +(UL*(SIGN*DRTA)-UL2*(SIGN*D1TTA)*GARO)*PPTA
     >         +(UL*(SIGN*D2TTA))*PPR
     >         +UL*(SIGN*D1TTA)*PPRTA
          P3RTA  = ((SIGN*D1TTA)*(-UL4*6*GARO*GARO*GARO
     >         +UL3*6*GARO*GA2RO -UL2*GA3RO) +3*(UL3*2*GARO*GARO
     >         - UL2*GA2RO)*(SIGN*DRTA) -UL2*3*GARO*(SIGN*D2RTA)
     >         +UL*(SIGN*D3RTA))*PP
     >         +3*((SIGN*D1TTA)*(UL3*2*GARO*GARO -UL2*GA2RO)
     >         -UL2*2*GARO*(SIGN*DRTA) +UL*(SIGN*D2RTA))*PPR
     >         +3*(-UL2*(SIGN*D1TTA)*GARO +UL*(SIGN*DRTA))*PP2R
     >         +UL*(SIGN*D1TTA)*PP3R
          PR3TA  = (UL*(SIGN*DR3TA)-UL2*GARO*(SIGN*D3TTA))*PP
     >         +(2*UL*(SIGN*DR2TA)-2*UL2*GARO*(SIGN*D2TTA))*PPTA
     >         +UL*(SIGN*D3TTA)*PPR
     >         +(UL*(SIGN*DRTA)-UL2*(SIGN*D1TTA)*GARO)*PP2TA
     >         +2*UL*(SIGN*D2TTA)*PPRTA
     >         +UL*(SIGN*D1TTA)*PPR2TA
          P2R2TA = ((SIGN*D2TTA)*(2*UL3*GARO*GARO-UL2*GA2RO)
     >         -2*UL2*GARO*(SIGN*DR2TA)+UL*(SIGN*D2R2TA))*PP
     >         +2*(-UL2*GARO*(SIGN*D2TTA)+UL*(SIGN*DR2TA))*PPR
     >         +((SIGN*D1TTA)*(2*UL3*GARO*GARO-UL2*GA2RO)
     >         -2*UL2*GARO*(SIGN*DRTA)+UL*(SIGN*D2RTA))*PPTA
     >         +(UL*(SIGN*D2TTA))*PP2R
     >         +2*(-UL2*(SIGN*D1TTA)*GARO+UL*(SIGN*DRTA))*PPRTA
     >         +UL*(SIGN*D1TTA)*PP2RTA

          IF(GAP .EQ. 0.D0) THEN
            F3RO   = 0.D0
            F4RO   = 0.D0
            F3TTA  = 0.D0
            F4TTA  = 0.D0
            F2RTA  = 0.D0
            FR2TA  = 0.D0
            F3RTA  = 0.D0
            FR3TA  = 0.D0
            F2R2TA = 0.D0
            IF(D.LE.0.D0) THEN
              F=1.D0
            ELSE
              F=0.D0
            ENDIF
          ELSE
            IF (P .GE. PLIM) THEN
              F3RO   = 0.D0
              F4RO   = 0.D0
              F3TTA  = 0.D0
              F4TTA  = 0.D0
              F2RTA  = 0.D0
              FR2TA  = 0.D0
              F3RTA  = 0.D0
              FR3TA  = 0.D0
              F2R2TA = 0.D0
            ELSEIF(P .LE. -PLIM) THEN
              F3RO   = 0.D0
              F4RO   = 0.D0
              F3TTA  = 0.D0
              F4TTA  = 0.D0
              F2RTA  = 0.D0
              FR2TA  = 0.D0
              F3RTA  = 0.D0
              FR3TA  = 0.D0
              F2R2TA = 0.D0
            ELSE
              EP2 = EP*EP
              EP3 = EP2*EP
              UEP  = 1.D0+EP
              UEP2 = UEP * UEP
              UEP3 = UEP2 * UEP
              UMEP = 1.D0-EP
              F4 = F3 * F
              F5 = F4 * F
              POL1 = (1.D0-4*EP+EP2)
              POL2 = (1.D0-10*EP+EP2)

              F3RO =  (EP* F4)*(-UEP2*P3RO-3*UEP*UMEP*P2RO*PRO
     >            -POL1*PRO2*PRO)
              F4RO = -(EP * F5)*(3*UEP2*UMEP*P2RO*P2RO+UEP3*P4RO
     >            +4*UMEP*UEP2*P3RO*PRO+6*UEP*POL1*P2RO*PRO2
     >            +UMEP*POL2*PRO2*PRO2)
              F3TTA = (EP* F4)*(-UEP2*P3TTA-3*UEP*UMEP*P2TTA*PTTA
     >            -POL1*PTTA2*PTTA)
              F4TTA = -(EP * F5)*(3*UEP2*UMEP*P2TTA*P2TTA+UEP3*P4TTA
     >            +4*UMEP*UEP2*P3TTA*PTTA+6*UEP*POL1*P2TTA*PTTA2
     >            +UMEP*POL2*PTTA2*PTTA2)
              F2RTA = -(EP * F4)*(UEP2*P2RTA+2*UMEP*UEP*PRO*PRTA
     >            +UMEP*UEP*P2RO*PTTA+POL1*PRO2*PTTA)
              FR2TA = -(EP * F4)*(UEP2*PR2TA+UMEP*UEP*P2TTA*PRO
     >            +2*UMEP*UEP*PRTA*PTTA+POL1*PRO*PTTA2)
              F3RTA = -(EP * F5)*(UEP3*P3RTA+3*UMEP*UEP2*(P2RTA*PRO
     >            +P2RO*PRTA)+3*UEP*POL1*PRO*(PRO*PRTA+P2RO*PTTA)
     >            +UMEP*UEP2*P3RO*PTTA+UMEP*POL2*PRO2*PRO*PTTA)
              FR3TA = -(EP * F5)*(UEP3*PR3TA+UMEP*UEP2*P3TTA*PRO
     >            +3*UMEP*UEP2*P2TTA*PRTA
     >            +3*UMEP*UEP2*PR2TA*PTTA
     >            +3*UEP*POL1*PTTA*(PRO*P2TTA+PRTA*PTTA)
     >            +UMEP*POL2*PRO*PTTA*PTTA2)
              F2R2TA =  -(EP * F5)*(UEP3*P2R2TA+UMEP*UEP2*P2RO*P2TTA
     >            +2*UMEP*UEP2*PR2TA*PRO+UEP*POL1
     >            *P2TTA*PRO2+2*UMEP*UEP2*(PRTA*PRTA+P2RTA*PTTA)
     >            +UEP*POL1*PTTA*(4*PRO*PRTA+P2RO*PTTA)
     >            +UMEP*POL2*PRO2*PTTA2)
            ENDIF
          ENDIF
       ENDIF

C************************ POSITION DE  M(X,Y) / FACE DE SORTIE **************
C     On se place dans la zone lineaire i.e IF( X.GE.CXS .AND. X.LE.AXS )

       XOS = SICOS * (Y - YBS) + XBS * SIN2S + X * COS2S
       YOS = SICOS * (X - XBS) + YBS * COS2S + Y * SIN2S
       YLS = YBS + ((X - XBS) * TANOS )
C       ... YL = COTE Y DE LA PROJECTION DE B SUR LA D1ROITE (MY)
C           (M=POINT COURANT), PARALLELEMENT A LA D1ROITE (ABC)
       DS = SQRT((X - XOS)**2 + (Y - YOS)**2 )
       DS2 = DS * DS
       ROOTS = DS
       IF(ROOTS.EQ.0.D0)
     > CALL ENDJOB('ERROR : roots=0 in dipsfa. Try DS=1e-99 in pgm',-99)
       SIGNS = 1.D0
       IF( Y .GE. YLS .OR. DS .LE. 1.D-6 ) THEN
         DS = -DS
         SIGNS = -1.D0
       ENDIF
C       DS=( DS + SHIFTE )

C*********************** CALCUL DES DERIVEES POUR LA FACE DE SORTIE ********
      GAPOS=LAMBDS(KMAG)

      C0 = CS(KMAG,1)
      C1 = CS(KMAG,2)
      C2 = CS(KMAG,3)
      C3 = CS(KMAG,4)
      C4 = CS(KMAG,5)
      C5 = CS(KMAG,6)

      XOROS    = SICOS*YRO     + COS2S*XRO
      XO2ROS   = SICOS*Y2RO    + COS2S*X2RO
      XOTTAS   = SICOS*YTTA    + COS2S*XTTA
      XO2TAS   = SICOS*Y2TTA   + COS2S*X2TTA
      XORTAS   = SICOS*YRTA    + COS2S*XRTA

      YOROS    = SICOS*XRO     + SIN2S*YRO
      YO2ROS   = SICOS*X2RO    + SIN2S*Y2RO
      YOTTAS   = SICOS*XTTA    + SIN2S*YTTA
      YO2TAS   = SICOS*X2TTA   + SIN2S*Y2TTA
      YORTAS   = SICOS*XRTA    + SIN2S*YRTA

      DXS    = X      -  XOS
      DXROS  = XRO    -  XOROS
      DX2ROS = X2RO   -  XO2ROS
      DXTTAS = XTTA   -  XOTTAS
      DX2TAS = X2TTA  -  XO2TAS

      DYS    = Y      -  YOS
      DYROS  = YRO    -  YOROS
      DY2ROS = Y2RO   -  YO2ROS
      DYTTAS = YTTA   -  YOTTAS
      DY2TAS = Y2TTA  -  YO2TAS

      DXRTAS   =  XRTA   -  XORTAS
      DYRTAS   =  YRTA   -  YORTAS

      ROOT2S  = ROOTS*ROOTS
      ROOT3S  = ROOT2S*ROOTS

      D1ROS   =   (DXROS*DXS  + DYROS*DYS)  / ROOTS
      D2ROS   = - (DXS*DXROS  + DYS*DYROS)**2 / ROOT3S
     >     + (DXROS*DXROS + DYROS*DYROS+DXS*DX2ROS+DYS*DY2ROS)
     >     / ROOTS
      D1TTAS  =   (DXTTAS*DXS + DYTTAS*DYS) / ROOTS
      D2TTAS  = - (DXS*DXTTAS + DYS*DYTTAS)**2 /ROOT3S
     >     + (DXTTAS*DXTTAS + DYTTAS*DYTTAS + DXS*DX2TAS
     >     + DYS*DY2TAS) / ROOTS
      DRTAS = - (DXS*DXTTAS + DYS*DYTTAS)*(DXS*DXROS + DYS*DYROS)
     >     /  ROOT3S
     >     + (DXTTAS*DXROS + DYTTAS*DYROS + DXS*DXRTAS +DYS*DYRTAS)
     >     /  ROOTS

      GAPS =  GAPOS * (RM/RO)**QGAPS
      GAPS2 = GAPS*GAPS
      GAPS3 = GAPS2*GAPS
      DSL = DS/GAPS

      ULS = 1.D0/GAPS
      ULS2 = ULS*ULS
      PS=C0+(C1+(C2+(C3+(C4+C5*DSL)*DSL)*DSL)*DSL)*DSL
      PP= C1+(2*C2+(3*C3+(4*C4+5*C5*DSL)*DSL)*DSL)*DSL
      PPP = 2*C2+(6*C3+(12*C4+20*C5*DSL)*DSL)*DSL
      EPS=EXP(PS)

      GAROS = -GAPOS*(QGAPS/RM)*(RM/RO)**(QGAPS+1)
      GA2ROS =  GAPOS*(QGAPS*(QGAPS+1.D0)/(RM*RM))
     >     *(RM/RO)**(QGAPS+2.D0)
      DDS =  ((SIGNS*D1ROS)*GAPS-DS*GAROS)/ GAPS2
      D2DS= ((SIGNS*D2ROS)*GAPS2-DS*GAPS*GA2ROS
     >     +2*GAROS*GAROS*DS-2*GAPS*(SIGNS*D1ROS)*GAROS)/GAPS3
      PPTAS  = ULS*(SIGNS*D1TTAS)*PPP
      PPRS   = DDS*PPP

      PROS    =  DDS * PP
      P2ROS   =  D2DS*PP + DDS*DDS*PPP
      PTTAS  = SIGNS * ULS*D1TTAS * PP
      P2TTAS = ULS*((SIGNS*D2TTAS)*PP + (SIGNS*D1TTAS)*PPTAS)
      PRTAS   = (-ULS2*GAROS*(SIGNS*D1TTAS)+ULS*(SIGNS*DRTAS))*PP
     >     +ULS*(SIGNS*D1TTAS)*PPRS

      IF(GAPS .EQ. 0.D0) THEN
        FROS    = 0.D0
        F2ROS   = 0.D0
        FTTAS   = 0.D0
        F2TTAS  = 0.D0
        FRTAS   = 0.D0
        IF(DS.LE.0.D0) THEN
          FS=1.D0
        ELSE
          FS=0.D0
        ENDIF
      ELSE
        IF  (PS .GE.  PLIM) THEN
          FS      = 0.D0
          FROS    = 0.D0
          F2ROS   = 0.D0
          FTTAS   = 0.D0
          F2TTAS  = 0.D0
          FRTAS   = 0.D0
        ELSEIF(PS .LE. -PLIM) THEN
          FS      = 1.D0
          FROS    = 0.D0
          F2ROS   = 0.D0
          FTTAS   = 0.D0
          F2TTAS  = 0.D0
          FRTAS   = 0.D0
        ELSE
          UEPS  = 1.D0+EPS
          UMEPS = 1.D0-EPS
          FS  = 1.D0/UEPS
          FS2 = FS * FS
          FS3 = FS2 * FS
          PROS2  = PROS*PROS
          PTTAS2 = PTTAS*PTTAS
          FROS =  -PROS*EPS * FS2
          F2ROS = -(EPS * FS3)*(P2ROS*UEPS+PROS2*UMEPS)
          FTTAS = -PTTAS*EPS * FS2
          F2TTAS= -(EPS * FS3)*(P2TTAS*UEPS+PTTAS**2*UMEPS)
          FRTAS  = -(EPS * FS3)*(PRTAS*UEPS+(PTTAS*PROS*UMEPS))
        ENDIF
      ENDIF

C      IF (KIRD.EQ.1) THEN
      IF (IDB.GE.4) THEN
         XO3ROS   = SICOS*Y3RO    + COS2S*X3RO
         XO4ROS   = SICOS*Y4RO    + COS2S*X4RO
         XO3TAS   = SICOS*Y3TTA   + COS2S*X3TTA
         XO4TAS   = SICOS*Y4TTA   + COS2S*X4TTA
         XO2RTAS  = SICOS*Y2RTA   + COS2S*X2RTA
         XOR2TAS  = SICOS*YR2TA   + COS2S*XR2TA
         XO3RTAS  = SICOS*Y3RTA   + COS2S*X3RTA
         XOR3TAS  = SICOS*YR3TA   + COS2S*XR3TA
         XO2R2TAS = SICOS*Y2R2TA  + COS2S*X2R2TA

         YO3ROS   = SICOS*X3RO    + SIN2S*Y3RO
         YO4ROS   = SICOS*X4RO    + SIN2S*Y4RO
         YO3TAS   = SICOS*X3TTA   + SIN2S*Y3TTA
         YO4TAS   = SICOS*X4TTA   + SIN2S*Y4TTA
         YO2RTAS  = SICOS*X2RTA   + SIN2S*Y2RTA
         YOR2TAS  = SICOS*XR2TA   + SIN2S*YR2TA
         YO3RTAS  = SICOS*X3RTA   + SIN2S*Y3RTA
         YOR3TAS  = SICOS*XR3TA   + SIN2S*YR3TA
         YO2R2TAS = SICOS*X2R2TA  + SIN2S*Y2R2TA

         DX3ROS = X3RO   -  XO3ROS
         DX4ROS = X4RO   -  XO4ROS
         DX3TAS = X3TTA  -  XO3TAS
         DX4TAS = X4TTA  -  XO4TAS

         DY3ROS = Y3RO   -  YO3ROS
         DY4ROS = Y4RO   -  YO4ROS
         DY3TAS = Y3TTA  -  YO3TAS
         DY4TAS = Y4TTA  -  YO4TAS

         DX2RTAS  =  X2RTA  -  XO2RTAS
         DXR2TAS  =  XR2TA  -  XOR2TAS
         DX3RTAS  =  X3RTA  -  XO3RTAS
         DXR3TAS  =  XR3TA  -  XOR3TAS
         DX2R2TAS =  X2R2TA -  XO2R2TAS

         DY2RTAS  =  Y2RTA  -  YO2RTAS
         DYR2TAS  =  YR2TA  -  YOR2TAS
         DY3RTAS  =  Y3RTA  -  YO3RTAS
         DYR3TAS  =  YR3TA  -  YOR3TAS
         DY2R2TAS =  Y2R2TA -  YO2R2TAS

         ROOT5S  = ROOT2S*ROOT3S
         ROOT7S  = ROOT5S*ROOT2S

         D3ROS   =   (3*(DXS*DXROS + DYS*DYROS)**3) / ROOT5S
     >        - (3*(DXS*DXROS + DYS*DYROS)*(DXROS**2 + DYROS**2
     >        + DXS*DX2ROS + DYS*DY2ROS)) / ROOT3S
     >        + (3*DXROS*DX2ROS + 3*DYROS*DY2ROS + DXS*DX3ROS
     >        + DYS*DY3ROS) / ROOTS
         D4ROS   = - 15*(DXS* DXROS + DYS* DYROS)**4 / ROOT7S
     >        + 9*(DXS*DXROS + DYS*DYROS)**2*(2*DXROS**2
     >        + 2*DYROS**2 +  2*DXS*DX2ROS + 2*DYS*DY2ROS) / ROOT5S
     >        - 3*(DXROS**2 + DYROS**2 + DXS*DX2ROS+ DYS*DY2ROS)**2
     >        /  ROOT3S
     >        - 4*(DXS*DXROS + DYS*DYROS)*(3*DXROS*DX2ROS
     >        + 3*DYROS*DY2ROS +  DXS*DX3ROS + DYS*DY3ROS) / ROOT3S
     >        + (3*DX2ROS**2 + 3*DY2ROS**2 + 4*DXROS*DX3ROS
     >        + 4*DYROS*DY3ROS + DXS*DX4ROS + DYS*DY4ROS) / ROOTS
         D3TTAS  =   (3*(DXS*DXTTAS + DYS*DYTTAS)**3) / ROOT5S
     >        - (3*(DXS*DXTTAS + DYS*DYTTAS)*(DXTTAS**2 + DYTTAS**2
     >        + DXS*DX2TAS + DYS*DY2TAS)) / ROOT3S
     >        + (3*DXTTAS*DX2TAS + 3*DYTTAS*DY2TAS
     >        + DXS*DX3TAS+DYS*DY3TAS) / ROOTS
         D4TTAS   =  -15*(DXS* DXTTAS + DYS* DYTTAS)**4 / ROOT7S
     >        + 9*(DXS*DXTTAS + DYS*DYTTAS)**2*(2*DXTTAS**2
     >        + 2*DYTTAS**2 + 2*DXS*DX2TAS + 2*DYS*DY2TAS) / ROOT5S
     >        - 3*(DXTTAS**2 + DYTTAS**2 + DXS*DX2TAS
     >        + DYS*DY2TAS)**2 /  ROOT3S
     >        - 4*(DXS*DXTTAS+DYS*DYTTAS)*(3*DXTTAS*DX2TAS
     >        + 3*DYTTAS*DY2TAS + DXS*DX3TAS + DYS*DY3TAS) / ROOT3S
     >        + (3*DX2TAS**2 + 3*DY2TAS**2 + 4*DXTTAS*DX3TAS
     >        + 4*DYTTAS*DY3TAS + DXS*DX4TAS + DYS*DY4TAS) / ROOTS
         D2RTAS =   3*(DXS*DXTTAS + DYS*DYTTAS)*(DXS*DXROS+DYS*DYROS)**2
     >        /  ROOT5S
     >        -  (DXS*DXROS + DYS*DYROS)*(2*DXTTAS*DXROS
     >        + 2*DYTTAS*DYROS + 2*DXS*DXRTAS + 2*DYS*DYRTAS) / ROOT3S
     >        -  ((DXS*DXTTAS + DYS*DYTTAS)*(DXROS**2 + DYROS**2
     >        + DXS*DX2ROS +  DYS*DY2ROS)) / ROOT3S
     >        +  (2*DXROS*DXRTAS + 2*DYROS*DYRTAS + DXTTAS*DX2ROS
     >        +  DYTTAS*DY2ROS +  DXS*DX2RTAS + DYS*DY2RTAS) / ROOTS
         DR2TAS = - (DXS*DXTTAS + DYS*DYTTAS)*(2*DXS*DXRTAS
     >        + 2*DXROS*DXTTAS+2*DYS*DYRTAS+2*DYROS*DYTTAS) / ROOT3S
     >        + (DXS*DXROS+DYS*DYROS)*(3*(DXS*DXTTAS + DYS*DYTTAS)**2
     >        / ROOT5S
     >        - (DXS*DX2TAS+DXTTAS**2+DYS*DY2TAS+DYTTAS**2) / ROOT3S)
     >        + (DXS*DXR2TAS + DX2TAS*DXROS + 2*DXRTAS*DXTTAS
     >        + DYS*DYR2TAS + DY2TAS*DYROS + 2*DYRTAS*DYTTAS) / ROOTS
         D3RTAS = - 15*(DXS*DXROS +DYS*DYROS)**3*(DXS*DXTTAS+DYS*DYTTAS)
     >        / ROOT7S
     >        - (DXS*DX3ROS + 3*DX2ROS*DXROS + DYS*DY3ROS
     >        + 3*DY2ROS*DYROS)* (DXS*DXTTAS + DYS*DYTTAS) / ROOT3S
     >        + 9*(DXS*DXROS + DYS*DYROS)*(DXS*DX2ROS + DXROS**2
     >        + DYS*DY2ROS + DYROS**2)*(DXS*DXTTAS+DYS*DYTTAS)/ ROOT5S
     >        - 3*(DXS*DXROS + DYS*DYROS)*(DXS*DX2RTAS+2*DXROS*DXRTAS
     >        +  DX2ROS*DXTTAS + DYS*DY2RTAS + 2*DYROS*DYRTAS
     >        + DY2ROS*DYTTAS) /  ROOT3S
     >        + (DXS*DX3RTAS + 3*DX2RTAS*DXROS +3*DX2ROS*DXRTAS
     >        + DX3ROS*DXTTAS + DYS*DY3RTAS + 3*DY2RTAS*DYROS
     >        + 3*DY2ROS*DYRTAS+DY3ROS*DYTTAS) /  ROOTS
     >        +  9*(DXS*DXROS + DYS*DYROS)**2*(DXS*DXRTAS+DXROS*DXTTAS
     >        +  DYS*DYRTAS + DYROS*DYTTAS) / ROOT5S
     >        -  3*(DXS*DX2ROS + DXROS**2 + DYS*DY2ROS
     >        + DYROS**2)*(DXS*DXRTAS + DXROS*DXTTAS + DYS*DYRTAS
     >        + DYROS*DYTTAS) / ROOT3S
         DR3TAS = (DXS*DXR3TAS + DX3TAS*DXROS + 3*DX2TAS*DXRTAS
     >        + 3*DXR2TAS*DXTTAS + DYS*DYR3TAS + DY3TAS*DYROS
     >        + 3*DY2TAS*DYRTAS+3*DYR2TAS*DYTTAS) / ROOTS
     >        - 3*(DXS*DXTTAS + DYS*DYTTAS)*(DXS*DXR2TAS +DX2TAS*DXROS
     >        + 2*DXRTAS*DXTTAS + DYS*DYR2TAS + DY2TAS*DYROS
     >        + 2*DYRTAS*DYTTAS) / ROOT3S
     >        + 3*(DXS*DXRTAS + DXROS*DXTTAS + DYS*DYRTAS
     >        + DYROS*DYTTAS)*((3*(DXS*DXTTAS + DYS*DYTTAS)**2)/ROOT5S
     >        - (DXS*DX2TAS + DXTTAS**2 + DYS*DY2TAS + DYTTAS**2)
     >        / ROOT3S)
     >        + (DXS*DXROS + DYS*DYROS)*((-15*(DXS*DXTTAS
     >        + DYS*DYTTAS)**3) / ROOT7S
     >        - (DXS*DX3TAS + 3*DX2TAS*DXTTAS + DYS*DY3TAS
     >        + 3*DY2TAS*DYTTAS) /  ROOT3S
     >        + (9*(DXS*DXTTAS + DYS*DYTTAS)*(DXS*DX2TAS+DXTTAS**2
     >        + DYS*DY2TAS + DYTTAS**2))/ROOT5S)
         D2R2TAS =  -(DXS*DXTTAS + DYS*DYTTAS)*(2*DXS*DX2RTAS
     >        + 4*DXROS*DXRTAS + 2*DX2ROS*DXTTAS + 2*DYS*DY2RTAS
     >        + 4*DYROS*DYRTAS + 2*DY2ROS*DYTTAS) / ROOT3S
         D2R2TAS = D2R2TAS
     >        + 3*(DXS*DXROS +DYS*DYROS)*(2*DXS*DXTTAS+2*DYS*DYTTAS)
     >        * (2*DXS*DXRTAS + 2*DXROS*DXTTAS + 2*DYS*DYRTAS
     >        + 2*DYROS*DYTTAS) / ROOT5S
         D2R2TAS = D2R2TAS
     >        - (2*(DXS*DXRTAS + DXROS*DXTTAS
     >        + DYS*DYRTAS+DYROS*DYTTAS)**2
     >        + 2*(DXS*DXROS + DYS*DYROS)*(DXS*DXR2TAS +DX2TAS*DXROS
     >        + 2*DXRTAS*DXTTAS + DYS*DYR2TAS + DY2TAS*DYROS
     >        + 2*DYRTAS*DYTTAS))/ROOT3S
         D2R2TAS = D2R2TAS
     >        + (DXS*DXROS+DYS*DYROS)**2
     >        *(-15*(DXS*DXTTAS+DYS*DYTTAS)**2 / ROOT7S
     >        + 3*(DXS*DX2TAS + DXTTAS**2 + DYS*DY2TAS
     >        + DYTTAS**2)/ROOT5S)
     >        + (DXS*DX2ROS + DXROS**2 + DYS*DY2ROS + DYROS**2)
     >        * (3*(DXS*DXTTAS + DYS*DYTTAS)**2/ROOT5S
     >        - (DXS*DX2TAS + DXTTAS**2 + DYS*DY2TAS + DYTTAS**2)
     >        /  ROOT3S)
     >        +  (DXS*DX2R2TAS +DX2ROS*DX2TAS+2*DXR2TAS*DXROS
     >        + 2*DXRTAS**2+2*DX2RTAS*DXTTAS+DYS*DY2R2TAS
     >        + DY2ROS*DY2TAS + 2*DYR2TAS*DYROS + 2*DYRTAS**2
     >        + 2*DY2RTAS*DYTTAS) /  ROOTS

         ULS3 = ULS2*ULS
         ULS4 = ULS3*ULS
         PPPP  = 6*C3+(24*C4+60*C5*DSL)*DSL
         PPPPP = 24*C4+120*C5*DSL
         GAPS4 = GAPS3*GAPS
         GAPS5 = GAPS4*GAPS
         QGAPS2 = QGAPS*QGAPS
         QGAPS3 = QGAPS2*QGAPS

         GA3ROS = -GAPOS*QGAPS*((2.D0+3*QGAPS+QGAPS2)/(RM*RM*RM))
     >        *(RM/RO)**(QGAPS+3.D0)
         GA4ROS =  GAPOS*QGAPS*((6.D0+11*QGAPS+6*QGAPS2+QGAPS3)
     >        /(RM*RM*RM*RM))*(RM/RO)**(QGAPS+4.D0)
         D3DS= (6*GAROS**2*(GAPS*(SIGNS*D1ROS)-DS*GAROS)
     >        -3*GAPS2*(GAROS*(SIGNS*D2ROS)+(SIGNS*D1ROS)*GA2ROS)
     >        +6*DS*GAPS*GAROS*GA2ROS
     >        +GAPS2*(GAPS*(SIGNS*D3ROS)-DS*GA3ROS))/GAPS4
         D4DS= (6*GA2ROS*GAPS2*(DS*GA2ROS-(SIGNS*D2ROS)*GAPS)
     >        -GAPS3*(DS*GA4ROS-GAPS*(SIGNS*D4ROS))
     >        +24*(SIGNS*D1ROS)*GA2ROS*GAPS2*GAROS
     >        + 8*DS*GA3ROS*GAPS2*GAROS
     >        -12*GAPS*GAROS*GAROS*(3*DS*GA2ROS-(SIGNS*D2ROS)*GAPS)
     >        -4*GAPS3*((SIGNS*D1ROS)*GA3ROS+(SIGNS*D3ROS)*GAROS)
     >        -24*GAROS*GAROS*GAROS*((SIGNS*D1ROS)*GAPS-DS*GAROS))/GAPS5
         PPPTAS = ULS*(SIGNS*D1TTAS)*PPPP
         PP2TAS = ULS*((SIGNS*D2TTAS)*PPP
     >        + (SIGNS*D1TTAS)*ULS*(SIGNS*D1TTAS)*PPPP)
         PP3TAS = ULS*((SIGNS*D3TTAS)*PPP
     >        + 3*(SIGNS*D2TTAS)*ULS*(SIGNS*D1TTAS)*PPPP
     >        + ULS2*(SIGNS*D1TTAS)*(SIGNS*D1TTAS)*(SIGNS*D1TTAS)*PPPPP)
         PPPRS  = DDS*PPPP
         PP2RS  = D2DS*PPP + DDS*DDS*PPPP
         PPP2RS = D2DS*PPPP+ DDS*DDS*PPPPP
         PP3RS  = D3DS*PPP + 3*D2DS*DDS*PPPP + DDS*DDS*DDS*PPPPP
         PPRTAS = -ULS2*GAROS*(SIGNS*D1TTAS)*PPP
     >        +ULS*(SIGNS*DRTAS)*PPP
     >        +ULS*(SIGNS*D1TTAS)*DDS*PPPP
         PPPRTAS = -ULS2*GAROS*(SIGNS*D1TTAS)*PPPP
     >        +ULS*(SIGNS*DRTAS)*PPPP
     >        +ULS*(SIGNS*D1TTAS)*DDS*PPPPP
         PP2RTAS = (-ULS2*2*(SIGNS*DRTAS)*GAROS+ULS*(SIGNS*D2RTAS)
     >        +(SIGNS*D1TTAS)*(2*ULS3*GAROS*GAROS-ULS2*GA2ROS))*PPP
     >        + 2*(ULS*(SIGNS*DRTAS)-ULS2*(SIGNS*D1TTAS)*GAROS)*PPPRS
     >        + ULS*(SIGNS*D1TTAS)*PPP2RS
         PPR2TAS  = (ULS*(SIGNS*DR2TAS)-ULS2*GAROS*(SIGNS*D2TTAS))*PPP
     >        +(ULS*(SIGNS*DRTAS)-ULS2*(SIGNS*D1TTAS)*GAROS)*PPPTAS
     >        +(ULS*(SIGNS*D2TTAS))*PPPRS
     >        +ULS*(SIGNS*D1TTAS)*PPPRTAS

         P3ROS   =  D3DS*PP + 3*D2DS*DDS*PPP + DDS*DDS*DDS*PPPP
         P4ROS   =  D4DS*PP + 4*D3DS*DDS*PPP + 3*D2DS*D2DS*PPP
     >        + 6*D2DS*DDS*DDS*PPPP + DDS*DDS*DDS*DDS*PPPPP
         P2RTAS  = (-ULS2*2*(SIGNS*DRTAS)*GAROS+ULS*(SIGNS*D2RTAS)
     >        +(SIGNS*D1TTAS)*(2*ULS3*GAROS*GAROS-ULS2*GA2ROS))*PP
     >        +2*(ULS*(SIGNS*DRTAS)-ULS2*(SIGNS*D1TTAS)*GAROS)*PPRS
     >        + ULS*(SIGNS*D1TTAS)*PP2RS
         PR2TAS  = (ULS*(SIGNS*DR2TAS)-ULS2*GAROS*(SIGNS*D2TTAS))*PP
     >        +(ULS*(SIGNS*DRTAS)-ULS2*(SIGNS*D1TTAS)*GAROS)*PPTAS
     >        +(ULS*(SIGNS*D2TTAS))*PPRS
     >        +ULS*(SIGNS*D1TTAS)*PPRTAS
         P3RTAS  = ((SIGNS*D1TTAS)*(-ULS4*6*GAROS*GAROS*GAROS
     >        +ULS3*6*GAROS*GA2ROS -ULS2*GA3ROS)+3*(ULS3*2*GAROS*GAROS
     >        - ULS2*GA2ROS)*(SIGNS*DRTAS)-ULS2*3*GAROS*(SIGNS*D2RTAS)
     >        +ULS*(SIGNS*D3RTAS))*PP
     >        +3*((SIGNS*D1TTAS)*(ULS3*2*GAROS*GAROS -ULS2*GA2ROS)
     >        -ULS2*2*GAROS*(SIGNS*DRTAS) +ULS*(SIGNS*D2RTAS))*PPRS
     >        +3*(-ULS2*(SIGNS*D1TTAS)*GAROS +ULS*(SIGNS*DRTAS))*PP2RS
     >        +ULS*(SIGNS*D1TTAS)*PP3RS
         PR3TAS  = (ULS*(SIGNS*DR3TAS)-ULS2*GAROS*(SIGNS*D3TTAS))*PP
     >        +(2*ULS*(SIGNS*DR2TAS)-2*ULS2*GAROS*(SIGNS*D2TTAS))*PPTAS
     >        +ULS*(SIGNS*D3TTAS)*PPRS
     >        +(ULS*(SIGNS*DRTAS)-ULS2*(SIGNS*D1TTAS)*GAROS)*PP2TAS
     >        +2*ULS*(SIGNS*D2TTAS)*PPRTAS+ULS*(SIGNS*D1TTAS)*PPR2TAS
         P2R2TAS = ((SIGNS*D2TTAS)*(2*ULS3*GAROS*GAROS-ULS2*GA2ROS)
     >        -2*ULS2*GAROS*(SIGNS*DR2TAS)+ULS*(SIGNS*D2R2TAS))*PP
     >        +2*(-ULS2*GAROS*(SIGNS*D2TTAS)+ULS*(SIGNS*DR2TAS))*PPRS
     >        +((SIGNS*D1TTAS)*(2*ULS3*GAROS*GAROS-ULS2*GA2ROS)
     >        -2*ULS2*GAROS*(SIGNS*DRTAS)+ULS*(SIGNS*D2RTAS))*PPTAS
     >        +(ULS*(SIGNS*D2TTAS))*PP2RS+2*(-ULS2*(SIGNS*D1TTAS)*GAROS
     >        +ULS*(SIGNS*DRTAS))*PPRTAS+ULS*(SIGNS*D1TTAS)*PP2RTAS
         P3TTAS = ULS*((SIGNS*D3TTAS)*PP + 2*(SIGNS*D2TTAS)*PPTAS
     >        + (SIGNS*D1TTAS)*PP2TAS)
         P4TTAS = ULS*((SIGNS*D4TTAS)*PP + 3*(SIGNS*D3TTAS)*PPTAS
     >        + 3*(SIGNS*D2TTAS)*PP2TAS + (SIGNS*D1TTAS)*PP3TAS)

         IF(GAPS .EQ. 0.D0) THEN
           F3ROS   = 0.D0
           F4ROS   = 0.D0
           F3TTAS  = 0.D0
           F3TTAS  = 0.D0
           F4TTAS  = 0.D0
           F2RTAS  = 0.D0
           FR2TAS  = 0.D0
           F3RTAS  = 0.D0
           FR3TAS  = 0.D0
           F2R2TAS = 0.D0
           IF(DS.LE.0.D0) THEN
             FS=1.D0
           ELSE
             FS=0.D0
           ENDIF
         ELSE
           IF (PS .GE.  PLIM) THEN
             F3ROS   = 0.D0
             F4ROS   = 0.D0
             F3TTAS  = 0.D0
             F3TTAS  = 0.D0
             F4TTAS  = 0.D0
             F2RTAS  = 0.D0
             FR2TAS  = 0.D0
             F3RTAS  = 0.D0
             FR3TAS  = 0.D0
             F2R2TAS = 0.D0
           ELSEIF(PS .LE. -PLIM) THEN
             F3ROS   = 0.D0
             F4ROS   = 0.D0
             F3TTAS  = 0.D0
             F4TTAS  = 0.D0
             F2RTAS  = 0.D0
             FR2TAS  = 0.D0
             F3RTAS  = 0.D0
             FR3TAS  = 0.D0
             F2R2TAS = 0.D0
           ELSE
             EPS2 = EPS*EPS
             EPS3 = EPS2*EPS
             UEPS  = 1.D0+EPS
             UEPS2 = UEPS * UEPS
             UEPS3 = UEPS2 * UEPS
             UMEPS = 1.D0-EPS
             FS4 = FS3 * FS
             FS5 = FS4 * FS
             POLS1 = (1.D0-4*EPS+EPS2)
             POLS2 = (1.D0-10*EPS+EPS2)

             F3ROS =  (EPS* FS4)*(-UEPS2*P3ROS-3*UEPS*UMEPS*P2ROS*PROS
     >           -POLS1*PROS2*PROS)
             F4ROS = -(EPS * FS5)*(3*UEPS2*UMEPS*P2ROS*P2ROS+UEPS3*P4ROS
     >           +4*UMEPS*UEPS2*P3ROS*PROS+6*UEPS*POLS1*P2ROS*PROS2
     >           +UMEPS*POLS2*PROS2*PROS2)
             F3TTAS = (EPS*FS4)*(-UEPS2*P3TTAS-3*UEPS*UMEPS*P2TTAS*PTTAS
     >           -POLS1*PTTAS2*PTTAS)
             F4TTAS = -(EPS*FS5)*(3*UEPS2*UMEPS*P2TTAS*P2TTAS
     >           +UEPS3*P4TTAS+4*UMEPS*UEPS2*P3TTAS*PTTAS
     >           +6*UEPS*POLS1*P2TTAS*PTTAS2+UMEPS*POLS2*PTTAS2*PTTAS2)
             F2RTAS = -(EPS * FS4)*(UEPS2*P2RTAS+2*UMEPS*UEPS*PROS*PRTAS
     >           +UMEPS*UEPS*P2ROS*PTTAS+POLS1*PROS2*PTTAS)
             FR2TAS = -(EPS * FS4)*(UEPS2*PR2TAS+UMEPS*UEPS*P2TTAS*PROS
     >           +2*UMEPS*UEPS*PRTAS*PTTAS+POLS1*PROS*PTTAS2)
             F3RTAS =-(EPS*FS5)*(UEPS3*P3RTAS+3*UMEPS*UEPS2*(P2RTAS*PROS
     >           +P2ROS*PRTAS)+3*UEPS*POLS1*PROS*(PROS*PRTAS
     >           +P2ROS*PTTAS)+UMEPS*UEPS2*P3ROS*PTTAS
     >           +UMEPS*POLS2*PROS2*PROS*PTTAS)
             FR3TAS = -(EPS * FS5)*(UEPS3*PR3TAS+UMEPS*UEPS2*P3TTAS*PROS
     >           +3*UMEPS*UEPS2*P2TTAS*PRTAS
     >           +3*UMEPS*UEPS2*PR2TAS*PTTAS
     >           +3*UEPS*POLS1*PTTAS*(PROS*P2TTAS+PRTAS*PTTAS)
     >           +UMEPS*POLS2*PROS*PTTAS*PTTAS2)
             F2R2TAS =-(EPS*FS5)*(UEPS3*P2R2TAS+UMEPS*UEPS2*P2ROS*P2TTAS
     >           +2*UMEPS*UEPS2*PR2TAS*PROS+UEPS*POLS1*P2TTAS*PROS2
     >           +2*UMEPS*UEPS2*(PRTAS*PRTAS+P2RTAS*PTTAS)
     >           +UEPS*POLS1*PTTAS*(4*PROS*PRTAS+P2ROS*PTTAS)
     >           +UMEPS*POLS2*PROS2*PTTAS2)
           ENDIF
         ENDIF
      ENDIF

C*************** CALCUL DE FL - Not used *************************

      FL      = 1.D0
      FROL    = 0.D0
      F2ROL   = 0.D0
      F3ROL   = 0.D0
      F4ROL   = 0.D0
      FTTAL   = 0.D0
      F2TTAL  = 0.D0
      F3TTAL  = 0.D0
      F4TTAL  = 0.D0
      FRTAL   = 0.D0
      F2RTAL  = 0.D0
      FR2TAL  = 0.D0
      F3RTAL  = 0.D0
      FR3TAL  = 0.D0
      F2R2TAL = 0.D0

C-------------------------------------------------------------
C    Calcul de B et des derivees de B par rapport a RO  et TTA

      ROI   = (RO-RM)/RM

      IF    (ITYPF.EQ.0) THEN
        HRC = HO*(1.D0+(COEF1+(COEF2+(COEF3
     >          +(COEF4+COEF5*ROI)*ROI)*ROI)*ROI)*ROI)
        DHRC  = (HO/RM)*(COEF1+(2.D0*COEF2+(3.D0*COEF3
     >          +(4.D0*COEF4+5.D0*COEF5*ROI)*ROI)*ROI)*ROI)
        D2HRC = (HO/(RM*RM))*(2.D0*COEF2+(6.D0*COEF3
     >          +(12.D0*COEF4+20.D0*COEF5*ROI)*ROI)*ROI)
      ELSEIF(ITYPF.EQ.1) THEN
        HRC = (HO+(COEF1+(COEF2+(COEF3
     >          +(COEF4+COEF5*ROI)*ROI)*ROI)*ROI)*ROI)
        DHRC  = (1.D0/RM)*(COEF1+(2.D0*COEF2+(3.D0*COEF3
     >          +(4.D0*COEF4+5.D0*COEF5*ROI)*ROI)*ROI)*ROI)
        D2HRC = (1.D0/(RM*RM))*(2.D0*COEF2+(6.D0*COEF3
     >          +(12.D0*COEF4+20.D0*COEF5*ROI)*ROI)*ROI)
      ENDIF

      R11=1.D0/RO
      R12=R11*R11

      FAC    =  F*FS*FL
      FACRO  =  FL*FS*FRO  +F*FS*FROL  + F*FL*FROS
      FAC2RO = 2.D0*(FL*FRO+F*FROL)*FROS + FS*(2.D0*FRO*FROL +
     >      FL*F2RO+F*F2ROL)    + F*FL*F2ROS
      FACTTA =  FL*FS*FTTA +F*FS*FTTAL + F*FL*FTTAS
      FAC2TA = 2.D0*(FL*FTTA+F*FTTAL)*FTTAS + FS*(2.D0*FTTA*FTTAL +
     >      FL*F2TTA    + F*F2TTAL)+F*FL*F2TTAS
      FACRTA  =  FRTA*FS*FL+FRO*FTTAS*FL+FRO*FS*FTTAL
     >          + FTTA*FROS*FL+F*FRTAS*FL+F*FROS*FTTAL
     >          + FTTA*FS*FROL+F*FTTAS*FROL+F*FS*FRTAL

      B      =  FAC*HRC
      BRO    =  FAC*DHRC  +   FACRO*HRC
      B2RO   =  FAC*D2HRC + 2.D0*FACRO*DHRC  +   FAC2RO*HRC
      BTTA   =  FACTTA*HRC
      B2TTA  =  FAC2TA*HRC
      BRTA   =  FACTTA*DHRC  +   FACRTA*HRC


c       iunit=88
c       if(ipass.GT.2000 .and. ipass.lt.2500)
c     > write(iunit,*) ' dipsi 1/qbr bri before :',1.d0/qbr,bri,ipass
           bri = 1.d0/qbr
c       if(ipass.GT.2000 .and. ipass.lt.2500)
c     >     write(iunit,*) ' dipsi bri after  :',1.d0/qbr,bri,ipass
c        call flush2(iunit,.false.)

      BZ     = B     * BRI
      BZX    = BTTA  * BRI
      BZY    = BRO   * BRI
      BZXX   = B2TTA * BRI
      BZYY   = B2RO  * BRI
      BZXY   = BRTA  * BRI

      BZX   =   BZX*R11
      BZXX  = ( BZXX*R11+BZY)*R11
      BZXY  = ( BZXY-BZX)*R11

C      IF (KIRD.EQ.1) THEN
      IF (IDB.GE.4) THEN

        IF    (ITYPF.EQ.0) THEN
         D3HRC =  (HO/(RM*RM*RM))*(6*COEF3+(24*COEF4+60*COEF5*ROI)*ROI)
         D4HRC =  (HO/(RM*RM*RM*RM))*(24*COEF4+120*COEF5*ROI)
        ELSEIF(ITYPF.EQ.1) THEN
         D3HRC =  (1.D0/(RM*RM*RM))*(6.D0*COEF3+
     >               (24.D0*COEF4+60.D0*COEF5*ROI)*ROI)
         D4HRC =  (1.D0/(RM*RM*RM*RM))*(24.D0*COEF4+120.D0*COEF5*ROI)
        ENDIF

         FAC3RO = 3*FROS*(2.D0*FRO*FROL+FL*F2RO+F*F2ROL)
     >        + 3*(FL*FRO+F*FROL)*F2ROS
     >        + FS*(3*FROL*F2RO+3*FRO*F2ROL+FL*F3RO+F*F3ROL)
     >        + F*FL*F3ROS
         FAC4RO =  6*(2.D0*FRO*FROL+FL*F2RO+F*F2ROL)*F2ROS
     >        + 4*FROS*(3*FROL*F2RO+3*FRO*F2ROL+FL*F3RO+F*F3ROL)
     >        + 4*(FL*FRO+F*FROL)*F3ROS
     >        + FS*(6*F2RO*F2ROL+4*FROL*F3RO+4*FRO*F3ROL+F4RO+F*F4ROL)
     >        + F*FL*F4ROS

         FAC3TA = 3*FTTAS*(2.D0*FTTA*FTTAL+FL*F2TTA+F*F2TTAL)
     >        + 3*(FL*FTTA+F*FTTAL)*F2TTAS
     >        + FS*(3*FTTAL*F2TTA+3*FTTA*F2TTAL+FL*F3TTA+F*F3TTAL)
     >        + F*FL*F3TTAS
         FAC4TA =  6*(2.D0*FTTA*FTTAL+FL*F2TTA+F*F2TTAL)*F2TTAS
C????     >        + 4*FTTAS*(3*FTTAL*F2TTA+3*FTTA*F2TTAL+FL*F3TTA F*F3TTAL)
     >        + 4*FTTAS*(3*FTTAL*F2TTA+3*FTTA*F2TTAL+FL*F3TTA+F*F3TTAL)
     >        + 4*(FL*FTTA+F*FTTAL)*F3TTAS
     >        + FS*(6*F2TTA*F2TTAL+4*FTTAL*F3TTA+4*FTTA*F3TTAL+FL*F4TTA
     >        + F*F4TTAL)
     >        + F*FL*F4TTAS
         FAC2RT =  2.D0*FROS*(FTTAL*FRO + FTTA*FROL + FL*FRTA + F*FRTAL)
     >       + 2.D0*(FL*FRO+F*FROL)*FRTAS + FTTAS*(2.D0*FRO*FROL+FL*F2RO
     >        +F*F2ROL)
     >        + FL*FTTA*F2ROS + F*FTTAL*F2ROS
     >        + FS*(2.D0*FROL*FRTA+2.D0*FRO*FRTAL+FTTAL*F2RO+FTTA*F2ROL
     >        + FL*F2RTA+F*F2RTAL)
     >        + F*FL*F2RTAS
         FACR2T = (2.D0*FTTAL*FTTAS+FS*F2TTAL+FL*F2TTAS)*FRO
     >        +(2.D0*FTTA*FTTAS+FS*F2TTA+F*F2TTAS)*FROL
     >        + (2.D0*FTTA*FTTAL+FL*F2TTA+F*F2TTAL)*FROS
     >        + 2.D0*(FS*FTTAL+FL*FTTAS)*FRTA
     >     + 2.D0*(FS*FTTA+F*FTTAS)*FRTAL + 2.D0*(FL*FTTA+F*FTTAL)*FRTAS
     >        + FL*FS*FR2TA + F*FS*FR2TAL + F*FL*FR2TAS
         FAC3RT = 3*FRTAS*(2.D0*FRO*FROL + FL*F2RO + F*F2ROL)
     >        + 3*(FTTAL*FRO+FTTA*FROL+FL*FRTA+F*FRTAL)*F2ROS
     >     + 3*FROS*(2.D0*FROL*FRTA+2.D0*FRO*FRTAL+FTTAL*F2RO+FTTA*F2ROL
     >        + FL*F2RTA + F*F2RTAL)
     >        + 3*(FL*FRO+F*FROL)*F2RTAS
     >        + FTTAS*(3*FROL*F2RO+3*FRO*F2ROL+FL*F3RO+F*F3ROL)
     >        + FL*FTTA*F3ROS + F*FTTAL*F3ROS
     >        + FS*(3*FRTAL*F2RO+3*FRTA*F2ROL+3*FROL*F2RTA
     >        + 3*FRO*F2RTAL+FTTAL*F3RO+FTTA*F3ROL+FL*F3RTA+F*F3RTAL)
     >        + F*FL*F3RTAS
         FACR3T=(3*FTTAS*F2TTAL+3*FTTAL*F2TTAS+FS*F3TTAL+FL*F3TTAS)*FRO
     >        +(3*FTTAS*F2TTA+3*FTTA*F2TTAS+FS*F3TTA+F*F3TTAS)*FROL
     >        +(3*FTTAL*F2TTA+3*FTTA*F2TTAL+FL*F3TTA+F*F3TTAL)*FROS
     >        + 3*(2.D0*FTTAL*FTTAS+FS*F2TTAL+FL*F2TTAS)*FRTA
     >        + 3*(2.D0*FTTA*FTTAS+FS*F2TTA+F*F2TTAS)*FRTAL
     >        + 3*(2.D0*FTTA*FTTAL+FL*F2TTA+F*F2TTAL)*FRTAS
     >        + 3*(FS*FTTAL+FL*FTTAS)*FR2TA
     >        + 3*(FS*FTTA+F*FTTAS)*FR2TAL
     >        + 3*(FL*FTTA+F*FTTAL)*FR2TAS + FL*FS*FR3TA + F*FS*FR3TAL
     >        + F*FL*FR3TAS
         FA2R2TA =  4*(FTTAL*FRO+FTTA*FROL+FL*FRTA+F*FRTAL)*FRTAS
     >        + 2.D0*FROS*(F2TTAL*FRO+F2TTA*FROL+2.D0*FTTAL*FRTA
     >        +2.D0*FTTA*FRTAL + FL*FR2TA + F*FR2TAL)
     >     + 2.D0*(FL*FRO+F*FROL)*FR2TAS + F2TTAS*(2.D0*FRO*FROL+FL*F2RO
     >        + F*F2ROL)
     >        + (2.D0*FTTA*FTTAL+FL*F2TTA+F*F2TTAL)*F2ROS
     > + 2.D0*FTTAS*(2.D0*FROL*FRTA+2.D0*FRO*FRTAL+FTTAL*F2RO+FTTA*F2ROL
     >        + FL*F2RTA + F*F2RTAL)
     >        + 2.D0*(FL*FTTA + F*FTTAL)*F2RTAS
     >    +FS*(4*FRTA*FRTAL+2.D0*FROL*FR2TA+2.D0*FRO*FR2TAL+ F2TTAL*F2RO
     >      +F2TTA*F2ROL+2.D0*FTTAL*F2RTA + 2.D0*FTTA*F2RTAL + FL*F2R2TA
     >        + F*F2R2TAL)+ F*FL*F2R2TAS

         B3RO = FAC*D3HRC + 3*FACRO*D2HRC+3*FAC2RO*DHRC  + FAC3RO*HRC
         B4RO = FAC*D4HRC + 4*FACRO*D3HRC+6*FAC2RO*D2HRC + 4*FAC3RO*DHRC
     >        + HRC*FAC4RO
         B3TTA  =  FAC3TA*HRC
         B4TTA  =  FAC4TA*HRC
         B2RTA  =  FACTTA*D2HRC + 2.D0*FACRTA*DHRC +  FAC2RT*HRC
         BR2TA  =  FAC2TA*DHRC  +   FACR2T*HRC
         B3RTA  =  FACTTA*D3HRC + 3*FACRTA*D2HRC+3*FAC2RT*DHRC
     >        + FAC3RT*HRC
         BR3TA  =  FAC3TA*DHRC  +   FACR3T*HRC
         B2R2TA =  FAC2TA*D2HRC + 2.D0*FACR2T*DHRC+  FA2R2TA*HRC

         BZXXX  = B3TTA  * BRI
         BZX4   = B4TTA  * BRI
         BZYYY  = B3RO   * BRI
         BZY4   = B4RO   * BRI
         BZXXY  = BR2TA  * BRI
         BZXYY  = B2RTA  * BRI
         BZX3Y  = BR3TA  * BRI
         BZXY3  = B3RTA  * BRI
         BZX2Y2 = B2R2TA * BRI

         BZX4  =(((BZX4 -8.D0*BZXX)*R11+6.D0*BZXXY-3.D0*BZY)*R11+
     >    3.D0*BZYY)*R12
         BZX3Y =(( 6.D0*BZX-3.D0*BZXXX*R11+BZX3Y-8.D0*BZXY)*R11+
     >    3.D0*BZXYY)*R12
         BZX2Y2=(((6.D0*BZXX*R11-4.D0*BZXXY+2.D0*BZY)*R11-
     >    2.D0*BZYY+BZX2Y2)*R11+BZYYY)*R11
         BZXY3 =((6.D0*(BZXY-BZX)*R11-3.D0*BZXYY)*R11+BZXY3)*R11
         BZXXX = ( BZXXX*R11 + 3.D0*BZXY - 2.D0*BZX )*R12
         BZXXY =(( BZXXY - 2.D0*BZXX*R11 - BZY )*R11 + BZYY)*R11
         BZXYY =   BZXYY*R11 + 2.D0*( BZX - BZXY )*R12

      ENDIF

      IF(KMAG.GT.1) THEN
        BZ0(1,1) = BZ0(1,1) + BZ
        BZ0(2,1) = BZ0(2,1) + BZX
        BZ0(1,2) = BZ0(1,2) + BZY
        BZ0(3,1) = BZ0(3,1) + BZXX
        BZ0(2,2) = BZ0(2,2) + BZXY
        BZ0(1,3) = BZ0(1,3) + BZYY
        BZ0(4,1) = BZ0(4,1) + BZXXX
        BZ0(3,2) = BZ0(3,2) + BZXXY
        BZ0(2,3) = BZ0(2,3) + BZXYY
        BZ0(1,4) = BZ0(1,4) + BZYYY
        BZ0(5,1) = BZ0(5,1) + BZX4
        BZ0(4,2) = BZ0(4,2) + BZX3Y
        BZ0(3,3) = BZ0(3,3) + BZX2Y2
        BZ0(2,4) = BZ0(2,4) + BZXY3
        BZ0(1,5) = BZ0(1,5) + BZY4
      ELSEIF(KMAG.EQ.1) THEN
         BZ0(1,1) = BZ
         BZ0(2,1) = BZX
         BZ0(1,2) = BZY
         BZ0(3,1) = BZXX
         BZ0(2,2) = BZXY
         BZ0(1,3) = BZYY
         BZ0(4,1) = BZXXX
         BZ0(3,2) = BZXXY
         BZ0(2,3) = BZXYY
         BZ0(1,4) = BZYYY
         BZ0(5,1) = BZX4
         BZ0(4,2) = BZX3Y
         BZ0(3,3) = BZX2Y2
         BZ0(2,4) = BZXY3
         BZ0(1,5) = BZY4
      ELSE
        STOP ' *** Error, SBR DIPSFA -> KMAG value'
      ENDIF

      IF (KMAG.LT.NBMAG) GOTO 30
      RETURN

 99   STOP
      END
