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
      SUBROUTINE CHAMC(X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     --------------------------------------------------------------
C     Champ CALCULE. APPELE PAR INTEGR A CHAQUE PAS D'INTEGRATION.
C     CALCULE, A PARTIR D'UNE FORMULE, LE Champ ET SES DERIVEES DANS
C     LE PLAN MEDIAN, AU POINT COURANT X,Y,Z=0.
C     --------------------------------------------------------------
      COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      PARAMETER(MCOEF=6)
      COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      COMMON/CHAMP/ BZ0(5,5), EZ0(5,5)
      COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QEL,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      COMMON/DDBXYZ/ DB(3,3),DDB(3,3,3)
      COMMON/D3BXYZ/ D3BX(3,3,3), D3BY(3,3,3), D3BZ(3,3,3)
      COMMON/D4BXYZ/ D4BX(3,3,3,3) ,D4BY(3,3,3,3) ,D4BZ(3,3,3,3)
      COMMON/DDEXYZ/ DE(3,3),DDE(3,3,3)
      COMMON/D3EXYZ/ D3EX(3,3,3), D3EY(3,3,3), D3EZ(3,3,3)
      COMMON/D4EXYZ/ D4EX(3,3,3,3) ,D4EY(3,3,3,3) ,D4EZ(3,3,3,3)
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      COMMON/MARK/ KART,KALC,KERC,KUASEX
      PARAMETER(MPOL=10)
      COMMON/MULTPE/ EM(MPOL),QLE(MPOL),QLS(MPOL),QFM(MPOL)
     >,QE(MPOL,6),QS(MPOL,6),RTQ(MPOL)
      COMMON/MULTPL/ BM(MPOL),DLE(MPOL),DLS(MPOL)
     >,DI(MPOL,10),DS(MPOL,10),RTB(MPOL)
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      COMMON/RIGID/ BORO,DPREF,DP,BR
 
      DIMENSION BT(5,15),ET(5,15)
      SAVE BT,ET

      DIMENSION FTAB(5,5)

      SAVE IRD2
 
CC     ** SI ON UTILISE CHAMBR ( => CALCUL ACCEPTANCE ) :
C      IF(LIMIT .EQ. 1) THEN
CC       ** FMAG > ( 1/(1+EXP(Co) )  INDIQUE QUE LA PARTICULE
CC          A FRANCHI LA FACE MAGNETIQUE .
C        IF(KALC .EQ. 3) THEN
C          IF(X .GE. XE .AND. X .LE. XS ) THEN
C            FMAG = 1.
C          ELSE
C            FMAG = 0.
C          ENDIF
C        ELSEIF(KALC .EQ. 1) THEN
C          FMAG = 1.
C        ENDIF
C      ENDIF
 
      IF( KALC .EQ. 3) THEN
C------- DIP, QUAD, SEXTU, ..., MULT, SOLENO...
 
        CALL RAZDRV(KFLD)
 
        IF(KUASEX.LE.MPOL+1) THEN
C--------- DIP, QUAD, SEXTU, etc. up to 20-POLE ; MULTIPOL ; ELMULT ; EBMULT

          IF(KFLD .EQ. MG) THEN
C----------- Magnetic
c        write(*,*) ' sbr chamc di ds ',di(1,1),di(1,2),di(1,3),di(1,4)
            CALL MULTIP(IDB,MG,KUASEX,X,Y,Z,BM,DLE,DLS,DI,DS,RTB,
     >               XE,XS,CE, CS,
     >               B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ,BT)
          ELSEIF(KFLD .EQ. LC) THEN
C----------- Electric
            CALL MULTIP(IDE,LC,KUASEX,X,Y,Z,EM,QLE,QLS,QE,QS,RTQ,
     >               XE,XS,QCE,QCS,
     >               E,DE,DDE,D3EX,D3EY,D3EZ,D4EX,D4EY,D4EZ,ET)
          ELSEIF(KFLD .EQ. ML) THEN
C----------- Electric & Magnetic
            CALL MULTIP(IDE,LC,KUASEX,X,Y,Z,EM,QLE,QLS,QE,QS,RTQ,
     >               XE,XS,QCE,QCS,
     >               E,DE,DDE,D3EX,D3EY,D3EZ,D4EX,D4EY,D4EZ,ET)
            CALL MULTIP(IDB,MG,KUASEX,X,Y,Z,BM,DLE,DLS,DI,DS,RTB,
     >               XE,XS, CE, CS,
     >               B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ,BT)
          ENDIF
       
        ELSEIF(KUASEX .EQ. 20) THEN
C--------- SOLENOID
          CALL SOLENF(X,Y,Z)
 
        ELSEIF(KUASEX .EQ. 21) THEN
C--------- WIENFILTER
C--------- Electric component

          Y0=Y
          Z0=Z
          CALL ROTX(EM(6),Y0,Z0)
C FM CeeRainer 11/2003. QCE and QCS where never used due to this bug, FF coeffs where always CE, CS. 
C          CALL BENDF(EM(1),MPOL,QLE,QLS,QE,QS,RTQ,X,Y0,
          CALL BENDF(EM(1),MPOL,XE,XS,QCE,QCS,QLE,QLS,QE,QS,RTQ,X,Y0,
     >                                                              EZ0)
          CALL SYMMED(Z0,IDZ,EZ0
     >                         ,E,DE,DDE,D3EX,D3EY,D3EZ,D4EX,D4EY)
          CALL XROTB(EM(6),E,DE,DDE,D3EX,D3EY,D3EZ,D4EX,D4EY,D4EZ)
C--------- Magnetic component

          Y0=Y
          Z0=Z
          CALL ROTX(BM(6),Y0,Z0)
C FM CeeRainer 11/2003. QCE and QCS where never used due to this bug, FF coeffs where always CE, CS. 
C          CALL BENDF(BM(1),MPOL,DLE,DLS,DI,DS,RTB,X,Y0,
          CALL BENDF(BM(1),MPOL,XE,XS,CE,CS,DLE,DLS,DI,DS,RTB,X,Y0,
     >                                                             BZ0)
          CALL SYMMED(Z0,IDZ,BZ0,
     >                          B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY)
          CALL XROTB(BM(6),B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
 
        ELSEIF(KUASEX .EQ. 22) THEN
C--------- EL2TUB. ELECTROSTATIQ 2-TUBE
          CALL EL2TU(X,Y,Z,BR)
 
        ELSEIF(KUASEX .EQ. 23) THEN
C--------- UNIPOT. ELECTROSTATIQ 3-TUBE
          CALL UNIPO(X,Y,Z,BR)

        ELSEIF(KUASEX.EQ.24) THEN
C--------- ELCYLDEF
          CALL ELCYL(MPOL,QLE,QLS,QE,QS,X,Y, 
     >                                                    E)
 
        ELSEIF(KUASEX.EQ.25) THEN
C--------- ELMIR
          Y0=Y
          Z0=Z
          IF(EM(6).NE.0.D0) CALL ROTX(EM(6),Y0,Z0)
          CALL ELMIRF(X,Z0,BR,
     >                           E,DE,DDE)
          IF(EM(6).NE.0.D0)
     >      CALL XROTB(EM(6),E,DE,DDE,D3EX,D3EY,D3EZ,D4EX,D4EY,D4EZ)

        ELSEIF(KUASEX .EQ. 26 ) THEN
C--------- ELCMIR 
          CALL ELCMIF(Y,Z,BR,
     >                         E,DE,DDE)
              
        ELSEIF(KUASEX .EQ. 28 ) THEN
C--------- HELIX
          CALL HELIXF(X,Y,Z,BR,RO,BO,
     >                                  B)
        ELSEIF(KUASEX .EQ. 29) THEN
C--------- COILS
          CALL COILSF(X,Y,Z)
 
        ELSEIF(KUASEX .EQ. 30) THEN
C--------- UNDULATOR
          CALL UNDULF(BM(1),X,Z,
     >                            B,DB,DDB)
        ELSE
          STOP ' SBR CHAMC :  No such field  installed !'
        ENDIF
 
      ELSEIF(KALC .EQ. 1) THEN
C------ Various fields B(X,Y,0), assuming median plane symmetry, 
C       extrapolation at Z performed from mid-plane field by Taylor exp.

        BN=BO/BR
 
        IF(KUASEX .EQ. 1) THEN
C         *** Champ QUADRUPOLAIRE A ETENDUE RADIALE DIPOLAIRE
C             = Q-POLE  SPECIAL  SPES2  .  GRADIENT=BO
          IF(ABS(Y).GT.30.D0) THEN
            IF(ABS(Y).GT.50.D0) THEN
              IF(ABS(Y).GT.80.D0) THEN
                BZ =0D0
                BZY=0D0
              ELSE
                BZ =BN*30-(Y-50)*BN
                BZY=-BN
              ENDIF
            ELSE
              BZ=30.D0*BN
              BZY=0D0
            ENDIF
          ELSE
            BZ =BN*Y
            BZY=BN
          ENDIF
 
        ELSE IF(KUASEX .EQ. 3 .OR. KUASEX .EQ. 4) THEN
C         ****  QUADISEX. Champ DIP+QUAD+SEX (KUASEX = 3)
C               SEXQUAD.            QUAD+SEX (         4)
          YO = Y + YCE
          IF    (YO.GE.0.D0) THEN
            EB = EB1
            EG = EG1
          ELSEIF(YO.LT.0.D0) THEN
            EB = EB2
            EG = EG2
          ENDIF
          IF    (KUASEX .EQ. 3) THEN
            U = 1D0
          ELSEIF(KUASEX .EQ. 4) THEN
            U = 0D0
          ENDIF
          BZ    = BN * (U + YO * (EN + YO * (EB + YO * EG)))
          BZY   = BN * (EN + YO * (2.D0 * EB + 3.D0 * YO * EG))
          BZYY  = BN * (2.D0 * EB + 6.D0 * YO * EG)
          BZYYY = BN * 6.D0 * EG
 
        ELSEIF(KUASEX .EQ. 5 ) THEN
C         ****  Champ CONSTANT B0 DANS UN AIMANT RECTANGULAIRE
C               DU TYPE  VENUS (LONG.MAGN.=1.4M , LARG=+/-.85M)
 
          IF(Y*Y .LE. RO*RO) THEN
            BZ    = BN
          ELSE
            BZ = 0D0
          ENDIF

C------------ Simulation of a short bend (RS in LHC ring, 06/1995)
C             Lorentzien within +/-X(BSEUIL) such that B(+/-X(BSEUIL))=0. 

             U = X / RO                         
             BZ    = (BO / ( 1.D0 + U * U ) - EB1)/BR

CC             Parabolic
C             U = (X - XLIM) / XLIM
C             BZ    = BN * ( 1.D0 - U * U )
CC             Sine
C             U = (X - XLIM)
C             BZ    = BN * cos(PI*U/(2.D0*XLIM)) 


        ELSE IF(KUASEX .EQ. 6 ) THEN
C         ****  Champ CONSTANT B0 DANS UN AIMANT CIRCULAIRE
C               DU TYPE  PS170 (RAYON DU BORD MAGN=RO= .7m EN PRINCIPE))
 
C         VERSION CIBLE SUR LA PERIPHERIE DE L'AIMANT
           XMR=X-RO
C         VERSION CIBLE AU CENTRE AIMANT :
C          XMR=X
          IF(XMR*XMR+Y*Y .LE. RO*RO ) THEN
            BZ    = BN
          ELSE
C           ** ON COMMET UNE PETITE ERREUR SUR SBDL EN ENTREE ET SORTIE
C              ( PLUS FAIBLE QUAND 'PAS' DIMINUE )
            BZ    = 0D0
          ENDIF
 
        ELSE IF(KUASEX .EQ. 7 ) THEN
C--------- Champ TOROIDAL DANS UN RECTANGLE
          CALL TOROID(X,Y)
 
        ELSEIF(KUASEX .EQ. 8 ) THEN

C--------- BEND
          IF(BM(6) .NE. 0) THEN
            Y0=Y
            Z0=Z
            CALL ROTX(BM(6),Y0,Z0)
          ENDIF
          CALL BENDF(BM(1),MPOL,XE,XS,CE,CS,DLE,DLS,DI,DS,RTB,X,Y,
     >                                                            BZ0)

        ELSEIF(KUASEX .EQ. 27) THEN
C--------- FFAG
C  Equivalence (X,ANGLE), (Y,RADIUS)
          IF    (IRD2.EQ.0) THEN 
C  Compute FFAG field and derivatives from analytical model
            CALL FFAGFA(IDB,X,Y,
     >                          BZ0)
          ELSEIF(IRD2.EQ.1) THEN
C  Compute FFAG field  and derivatives from flying field-mesh
            CALL FFAGF(X,Y,
     >                     DA,DR,FTAB)
            AAA = ZERO
            RRR = ZERO
            CALL INTPLF(Y,AAA,RRR,DA,DR,FTAB,IRD, 
     >                                             BZ0)
          ELSE
            STOP '*** SBR CHAMC. No such field case'
          ENDIF


        ELSEIF(KUASEX .EQ. 30 ) THEN
C--------- Free!!

        ELSEIF(KUASEX .EQ. 31) THEN
C--------- DIPOLE
C  Equivalence (X,ANGLE), (Y,RADIUS)
          CALL DIPF(X,Y,
     >                  AAA,RRR,DA,DR,FTAB)
          CALL INTPLF(Y,AAA,RRR,DA,DR,FTAB,IRD,
     >                                           BZ0)

        ELSEIF(KUASEX .EQ. 32) THEN
C--------- DIPOLES
C  Equivalence (X,ANGLE), (Y,RADIUS)
          IF    (IRD2.EQ.0) THEN 
              CALL DIPSFA(IDB,X,Y,
     >                            BZ0)
          ELSEIF(IRD2.EQ.1) THEN
            CALL DIPSF(X,Y,
     >                     DA,DR,FTAB)
C     >                     AAA,RRR,DA,DR,FTAB)
            CALL INTPLF(Y,AAA,RRR,DA,DR,FTAB,IRD,
     >                                             BZ0)
          ENDIF

        ELSEIF(KUASEX .EQ. 33) THEN
C--------- FFAG-SPI
C  Equivalence (X,ANGLE), (Y,RADIUS)
          IF    (IRD2.EQ.0) THEN 
C  Compute FFAG field and derivatives from analytical model
            CALL FFGSPA(IDB,X,Y,
     >                          BZ0)
          ELSEIF(IRD2.EQ.1) THEN
C  Compute FFAG field  and derivatives from flying field-mesh
            CALL FFGSPF(X,Y,
     >                      DA,DR,FTAB)
            AAA = ZERO
            RRR = ZERO
            CALL INTPLF(Y,AAA,RRR,DA,DR,FTAB,IRD, 
     >                                             BZ0)
          ELSE
            STOP '*** SBR CHAMC. No such field case' 
          ENDIF

        ENDIF
C---------- ENDIF KUASEX
 
        CALL SYMMED(Z,IDZ,BZ0, 
     >                       B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY)

        IF(KUASEX .EQ. 8 ) THEN
          IF(BM(6) .NE. 0) THEN
            CALL XROTB(BM(6),B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
          ENDIF
        ENDIF

      ENDIF
C---------- ENDIF KALC
 
      IF    (KFLD .EQ. MG) THEN
        CALL DBDXYZ(IDB,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
      ELSEIF(KFLD .EQ. LC) THEN
        CALL DBDXYZ(IDE,DE,DDE,D3EX,D3EY,D3EZ,D4EX,D4EY,D4EZ)
      ELSEIF(KFLD .EQ. ML) THEN
        CALL DBDXYZ(IDE,DE,DDE,D3EX,D3EY,D3EZ,D4EX,D4EY,D4EZ)
        CALL DBDXYZ(IDB,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
      ENDIF
      RETURN

      ENTRY CHAMC6(IRD2I)
      IRD2 = IRD2I
      RETURN
      END
