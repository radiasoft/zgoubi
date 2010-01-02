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
      SUBROUTINE DIPCI(SCAL, 
     >                      DSREF,IRD,IDB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
C      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DPPP,QBR,BRI
  
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
      DIMENSION DXE(NM),THETA(NM),R1(NM),U1(NM),U2(NM),R2(NM) 
      DIMENSION DXS(NM),THETAS(NM),R1S(NM),U1S(NM),U2S(NM),R2S(NM) 
      DIMENSION DX3(NM),THETA3(NM),R13(NM),U13(NM),U23(NM),R23(NM) 

      SAVE ACN,DRM,HNORM,CE,CS,CC
      SAVE IND,CIND
      SAVE LAMBDE,QSIE,NCOEFE,SHIFTE
      SAVE LAMBDS,QSIS,NCOEFS,SHIFTS
      SAVE LAMBD3,QSI3,NCOEF3,SHIFT3
      SAVE DXE,THETA,R1,U1,U2,R2
      SAVE DXS,THETAS,R1S,U1S,U2S,R2S
      SAVE NN, RESOL, ITYPF
      SAVE NBMAG, NBFACE

      DIMENSION BZ0(5,5)

      CHARACTER TYPCAL(2)*14

      SAVE TYPCAL

      PARAMETER (PLIM=80.D0)
      LOGICAL SHARPE, SHARPS

      SAVE FINTE, FINTS

      DATA TYPCAL / ' analytic', ' interpolation'/
      DATA FINTE, FINTS / 0.2D0, 0.2D0 /

C  NBMAG=number of magnets.  AT=total extent angle of field 
      NP = 2
      NBMAG = NINT(A(NOEL,NP))
      NP=NP+1 
      AT    = A(NOEL,NP)
      NP=NP+1 
      RM    = A(NOEL,NP)

      IF(NRES.GT.0) WRITE(NRES,*) ' WARNING : SBR DIPC IS UNDER ',
     >  ' DEVELOPEMENT !!'

      IF(NRES.GT.0) WRITE(NRES,200) NBMAG,AT,RM
  200 FORMAT(20X,1P,'Dipole  magnet  N-tuple,  # of  units  N = ',I2,//,
     1 11X,' Total  extent  of  magnet  XT = ',G12.4,' degs.',
     2 /,11X,' Reference  Y value  of  sector  is  YM = ',G12.4,' cm')

C  HNORM=Champ MAX DANS LE DIPOLE.
C  CIND= field indices
C  AT=etendue totale de la zone de champ, XCENT='X AU CENTRE',
C    YM = Y  MOYEN DE LA zone DE Champ. 
C  NBFACE(KMAG)=(2)3 : dipole limited by (2)3 field boundaries

C-----  a list of NBMAG.LE.NM magnets
      KMAG = 0
 10   CONTINUE
      KMAG = KMAG+1
      NP=NP+1 
      XCENT = A(NOEL,NP)
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

      ACN(KMAG) = XCENT     
      IF(NRES.GT.0) THEN
        WRITE(NRES,100) KMAG,ACN(KMAG),DRM(KMAG),HNORM(KMAG)
 100    FORMAT(/,5X,'Dipole  magnet # ',I1,/,1P,
     >  11X,' Positionning  abscissa  XCENT : ',G12.4,' cm',/,
     >  11X,' Positionning  wrt.  YM  : ',G12.4,' cm',/,
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
      DXE(KMAG) = A(NOEL,NP)
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
     >  1P,G12.4,' cm ,  gap-k =',G12.4)
        WRITE(NRES,127) NCOEFE(KMAG),(CE(KMAG,I),I=1,6),SHIFTE(KMAG)
127     FORMAT (10X,' COEFFICIENTS :',I3,6F10.5
     2  ,/,10X,' Shift  of  EFB  = ',G12.4,' CM',/)
        IF(LAMBDE(KMAG) .LE. 0.D0)  WRITE(NRES,FMT='(/,
     >      ''  ***  Warning : sharp edge '',
     >      ''model entails vertical wedge focusing approximated with'',
     >        '' first order kick  ***'')')
        WRITE(NRES,103) DXE(KMAG),THETA(KMAG),R1(KMAG),U1(KMAG),
     >          U2(KMAG),R2(KMAG)
 103    FORMAT(10X,'DX_face =',F7.2,' cm',
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
      DXS(KMAG)  = A(NOEL,NP)
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
        IF(LAMBDS(KMAG) .LE. 0.D0)  WRITE(NRES,FMT='(/,
     >      ''  ***  Warning : sharp edge '',
     >      ''model entails vertical wedge focusing approximated with'',
     >        '' first order kick  ***'')')
        WRITE(NRES,103) DXS(KMAG),THETAS(KMAG),R1S(KMAG),U1S(KMAG),
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
      DX3(KMAG) = A(NOEL,NP)
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
113       FORMAT(20X,' Face centred on direction XCENT+OMEGA,  at ',
     >    G12.4,' CM'/)
          WRITE(NRES,103) DX3(KMAG),THETA3(KMAG),R13(KMAG),U13(KMAG),
     >         U23(KMAG),R23(KMAG)
          IF(R13(KMAG)*R23(KMAG) .EQ. 0.D0) WRITE(NRES,123)
        ENDIF
      ENDIF

      IF(KMAG.LT.NBMAG) GOTO 10
C-----------------------------
      
      SHARPE = .TRUE.
      SHARPS = .TRUE.
      DO 57 KMAG = 1, NBMAG
        SHARPE = SHARPE .AND. (LAMBDE(KMAG).LE.0.D0)
        SHARPS = SHARPS .AND. (LAMBDS(KMAG).LE.0.D0)
 57   CONTINUE
C------- Correction for wedge
      IF(NRES.GT.0) WRITE(NRES,*) ' WARNING, SBR DIPCI : ',
     >' Make sure you want hard edge correction with ',
     >' FINT at entrance / exit = ', QSIE(1), QSIS(1)
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
      IRDA = NINT(A(NOEL,NP))
C Get resol, or idb
      NP=NP+1 
      RESOL=A(NOEL,NP)
      IF    (IRDA.NE.0) THEN
C        interpolation method
        IF(SHARPE .OR. SHARPS) CALL ENDJOB
     >    ('ERROR :  sharp edge not compatible with num. deriv.',-99)
        IRD = IRDA
        IRDA = 1
        IF    (IRD.EQ.2) THEN 
          NN=3
        ELSEIF(IRD.EQ.25) THEN
          NN=5
        ELSEIF(IRD.EQ.4) THEN
          NN=5
        ELSE
          STOP ' *** ERROR - SBR DIP2, WRONG VALUE IRD'
        ENDIF
      ELSEIF(IRDA.EQ.0) THEN
C        analytic
C        IDB = NINT(10*A(NOEL,NP))
        IDB = NINT(RESOL)
        IF(IDB.NE.4) IDB=2
      ENDIF
      CALL CHAMC6(IRDA)

C Get flag field type. iord.option has the form xx.yy, iord=2, 25 or 4 and yy=01-99
C      ITYPF = NINT(100.D0* (A(NOEL,NP-1)-DBLE(IRD)) )
      ITYPF = NINT(100.D0* (A(NOEL,NP-1)-DBLE(NINT(A(NOEL,NP-1)))) )
      IF(NRES .GT. 0) THEN
        WRITE(NRES,FMT='(/,5X,'' itypf= '',i6)') itypf
        IF    (ITYPF.EQ.0) THEN
          WRITE(NRES,FMT='(/,5X,'' Field formula : default  '')')
        ELSEIF(ITYPF.EQ.1) THEN
          WRITE(NRES,FMT='(/,5X,'' Field formula : type 2  '')')
        ENDIF
      ENDIF

      AE=0.D0
      AS=0.D0
      AT = AT
      XI = 0.D0
      XF = AT

C Formula to be revisited...
      KMAG = 1
      DSREF = RM * (  (DXE(1)-DXS(KMAG)) + 
     >         TAN(ACN(1) - DXE(1)) + 
     >         TAN(AT - ACN(KMAG) + DXS(KMAG)) )

      IF(NRES.GT.0) THEN
        WRITE(NRES,FMT='(/,5X,'' Field & deriv. calculation :'',A)') 
     >  TYPCAL(IRDA+1)
        IF    (IRDA.NE.0) THEN
          IF(IRD .EQ. 2) THEN
            WRITE(NRES,121) '3*3', RESOL
 121        FORMAT(20X,A3,
     >        '-point  interpolation, size of flying mesh :  ',
     >        'integration step /',1P,G12.4)
          ELSEIF(IRD .GE. 4) THEN
C----------- IRD= 4 OR 25
            WRITE(NRES,121) '5*5', RESOL
          ENDIF
        ELSEIF(IRDA.EQ.0) THEN
          WRITE(NRES,FMT='(20X,
     >    ''Derivatives computed to order '',I1)') IDB
        ENDIF
      ENDIF

      RETURN

C--------------------------------------------------------------------
C----- Numerical calculation of B & derivatives ---------------------
      ENTRY DIPCF(TTA,RO,
     >                   DXG,DYG,FTAB)
      CALL INTEG5(
     >            STEP)
                  
      CALL RAZ(FTAB,5*5)
      DYG = STEP/RESOL
      DXG = DYG

      KMAG = 0
 20   CONTINUE
      KMAG = KMAG+1

      RRM(KMAG) = RM +DRM(KMAG)

      XMEG = DXE(KMAG)
      TETA = THETA(KMAG) * RAD
      UT = XMEG - TETA
      SINO = SIN(UT)
      COSO = COS(UT)
      TANO = SINO / COSO
 
C  CALCUL DES PARAMETRES DE FACE  :
C  AXE X DU REFERENTIEL  = PARALLELE  DIRCTN XCENT, SENS DES RAYONS CROISSANTS
C  AXE Y DU REFERENTIEL  = ORTHOGONAL DIRCTN XCENT, DIRIGE VERS FACE D'ENTREE
C  ORIGINE DU REFERENTIEL= A L'INTERSECTION DE RM ET XCENT
 
C  PROJECTIONS DU RAYON MOYEN SUR LES AXES X,Y
      XB = RRM(KMAG) * ( COS(XMEG) - 1.D0)
      YB = RRM(KMAG) * SIN(XMEG)
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
 
      XMEGS = DXS(KMAG)
      TETAS = THETAS(KMAG) * RAD
      UTS= XMEGS- TETAS
      SINOS= SIN(UTS)
      COSOS= COS(UTS)
      TANOS= SINOS/ COSOS
      XBS = RRM(KMAG) * ( COS(XMEGS) - 1.D0)
      YBS = RRM(KMAG) * SIN(XMEGS)
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
        XMEG3 = DX3(KMAG)
        TETA3 = THETA3(KMAG) * RAD
        UT3= XMEG3- TETA3
        IF(ABS(UT3) .GE. .5d0*PI) THEN
          IF(NRES.GT.0) WRITE(NRES,139)
139       FORMAT(/,20X,' EFB lateral, ATTENTION : ABS(UT3) > 90 deg =>',
     >    '  ambiguous  for  calculation  of  field')
        ENDIF
 
        SINO3= SIN(UT3)
        COSO3= COS(UT3)
        TANO3= SINO3/ COSO3
        XB3 = RM3(KMAG)* COS(XMEG3) - RM
        YB3 = RM3(KMAG)* SIN(XMEG3)
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
        ROJ = RO + DYG * DBLE(NN-JRO-INT(NN/2))

        DO 1 ITTA = 1,NN
          TTAI = TTA - DXG * DBLE(NN-ITTA-INT(NN/2))
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
 104        FORMAT(/,5X,10('*'),' ERREUR PARAMETRES FACE ENTREE',/)
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
114         FORMAT(/,5X,10('*'),' ERREUR PARAMETRES FACE SORTIE',/)
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
124             FORMAT(/,5X,9('*'),'ERREUR PARAMETRES FACE LATERALE',/)
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
 
C          ROI= (ROJ-RM)/RM
C          ROI= ROJ-RM
          ROI= ROJ-RRM(KMAG)

          IF    (ITYPF.EQ.0) THEN
C            ROI= ROI/RM
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
            SF = HNORM(KMAG)
            IF(IND(KMAG).NE.0) THEN
              PROI = 1.D-1 ! convert CIND from T/m^* to kG/cm^*
              DO KND = 1, IND(KMAG)
                PROI = PROI * ROI
                SF = SF + CIND(KMAG,KND) * PROI
              ENDDO
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
      ENTRY DIPCFA(IDB,TTA,RO,
     >                        BZ0)
      CALL ENDJOB('** SBR DIPCI : analytical model not implemented',-99)
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

      RETURN

 99   STOP
      END 
