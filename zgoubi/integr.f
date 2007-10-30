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
      SUBROUTINE INTEGR(EVNT,BACKW,MIRROR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EVNT,BACKW,MIRROR
C     -------------------------------------------------------------
C     APPELE PAR TRANSF.
C     CALCULE UNE TRAJECTOIRE DE X A XLIM; DXI=PAS EN X OU EN ANGLE
C     PAS = PAS EN S : cf. DEPLA, COFIN
C     -------------------------------------------------------------
      COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
      COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      COMMON/CONST/ CL9,CEL,PI,RAD,DEG,QE ,AMPROT, CM2M
C rustine RACCAM
C      INCLUDE 'MXLD.H'
C      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      COMMON/EFBS/ AFB(2), BFB(2), CFB(2), IFB
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      COMMON/RIGID/ BORO,DPREF,DP,BR
      COMMON/STEP/ TPAS(3), KPAS
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
 
      LOGICAL WEDGE, WEDGS
      SAVE WEDGE, WEDGS

      PARAMETER (DLIM = 1.D-6)
      INCLUDE 'MXSTEP.H'
      PARAMETER (PI24= 2.467401101D0)

      LOGICAL BACK, CHREG

      SAVE NSTEP
      DIMENSION  AREG(2),BREG(2),CREG(2), AREGI(2),BREGI(2),CREGI(2)
      SAVE AREG,BREG,CREG,KREG

      DATA WEDGE, WEDGS / .FALSE.,  .FALSE./
      DATA WDGE, WDGS, FFXTE, FFXTS / 0.D0, 0.D0, 0.D0, 0.D0 /
      SAVE WDGE, WDGS, FFXTE, FFXTS

C rustine RACCAM
c      logical first
c      data first / .true. /

      KDR = 2
      ISORT=1
      CT=COS(T)
      ST=SIN(T)

      IF(KPAS .NE. 0)  THEN
C------- 3-region step &/or variable step
        IF(TPAS(1) .EQ. 0.D0) THEN
C--------- Entrance is sharp edge
          KREG = 2
        ELSE
          KREG = 1
        ENDIF
        PAS=TPAS(KREG)
        DXI = PAS
      ENDIF

      NSTEP = 0
      X2 = X
      Y2 = Y  
      Z2 = Z

C----- DEBUT DE BOUCLE SUR DXI
    1 CONTINUE

C---------------------  Some tests to possibly stop integration 
      NSTEP = NSTEP+1
      IF(NSTEP .GT. MXSTEP) THEN
        WRITE(ABS(NRES),*) 'Maximum # steps allowed = ',MXSTEP
        WRITE(6,*) '      Maximum # steps allowed = ',MXSTEP
        CALL KSTOP(2,IT,KEX,*97)
      ENDIF

      BACK = T*T .GT. PI24
      IF(.NOT.BACKW) THEN
C------- PARTICLES BENT BACKWARD THOUGH NOT ALLOWED ARE STOPPED
        IF(BACK) CALL KSTOP(3,IT,KEX,*97)  
      ENDIF
C-------------------------------------------------------------------
 
      IF(KALC .EQ. 2 ) THEN
C-------  Compute  B(X,Y,0)  from  field  maps

        CALL CHAMK(X2,Y2,Z2,*99)

        IF(KERK .NE. 0) THEN
C          IF(ABS(BZ).GT.0.020)  KEX =-1
          IF(ISORT .EQ. 2)  GOTO 7
          IF(NRES.GT.0) WRITE(NRES,100)Y,X,T,IT,' HORS LIMITES'
  100     FORMAT (/,4X,'Y =',F8.2,' X =',F8.2,' T =',F6.3
     >    ,' TRAJECTOIRE ',I3,A,' Champ ')
          ISORT=2
          GOTO 7
        ENDIF
        IF(ISORT .EQ. 1)  GOTO 7
        IF(NRES.GT.0) WRITE(NRES,100) Y,X,T,IT,' REVENUE DANS'
        ISORT=1
 
      ELSE
C-------  Compute  B(X,Y,Z), E(X,Y,Z)  from  mathematical  2D or 3D  field  models

        CALL CHAMC(X2,Y2,Z2)

      ENDIF

 7    CONTINUE

      IF(NSTEP .EQ. 1) THEN
C-------- Entrance wedge correction in BEND, in MULTIPOL(if non zero B1), 
        IF(WEDGE) CALL WEDGKI(1,T,Z,P,WGDE,FFXTE)
      ENDIF

      CALL DEVTRA

      IF(KPAS.NE.0 ) THEN
C------- Test case diffrnt step values @ entrance|body|exit. 
C         Exist entrance and/or exit fringe field regions and XPAS is coded. 

C-------- CHREG is .true. if particle is going to next region 
        IF(CHREG(KART,X,Y,XE,AREG,BREG,CREG,DXI,TPAS,
     >                                               KREG,AL,BL,D)) THEN
            CL = -D
            ST = SIN(T)
            CT = COS(T)
            COSTA = AL*CT + BL*ST
            PAF = D/COSTA
            CALL DEPLA(PAF)
            COSTA=AL*CT+BL*ST
            CALL ITER(AL,BL,CL,PAF,COSTA,KEX,*97)
            CALL COFIN
     >        (KART,NPLT,LST,PAF,Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT,EVNT,
     >         *97)
            CT=COS(T)
            ST=SIN(T)
            PAS = TPAS(KREG)
            DXI = PAS

            GOTO 1
        ENDIF
      ENDIF

C----------------   TEST  DROITE(S) DE COUPURE SORTIE
C                   Used either in maps, 
C                   or in case of sharp edge exit in BEND, WIENF, UNDUL, ELMIR...
C FM 08/99      IF(IDRT .GE. 2) THEN
      IF(IDRT .GE. 1) THEN
C------- DROITE(S) DE COUPURE EN SORTIE

        IF( KART .EQ. 1) THEN
C--------- COORDONNEES CARTESIENNES
          D = ABS(AM(KDR)*X + BM(KDR)*Y + CM(KDR))

          IF(D .LT. DXI) THEN
            AL = AM(KDR)
            BL = BM(KDR)
            CL = -D
            ST = SIN(T)
            CT = COS(T)
            COSTA = AL*CT + BL*ST
            PAF = D/COSTA

            IF(.NOT.MIRROR) THEN 
C------------- Step onto the exit 'DROITE'

              CALL DEPLA(PAF)
              COSTA=AL*CT+BL*ST
              CALL ITER(AL,BL,CL,PAF,COSTA,KEX,*97)
C              CALL ITER(AL,BL,CL,PAF,COSTA,KEX,*99)
              CALL COFIN
     >          (KART,NPLT,LST,PAF,Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT,EVNT,
     >           *97)
C     >           *99)
              IF(LST .EQ. 2) THEN 
                CALL IMPPLA(NPLT,PAF,Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT)
              ELSEIF(LST .EQ. 3) THEN 
                CALL IMPPLA(  30,PAF,Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT)
              ENDIF

C FM 08/99              IF    (KDR .EQ. IDRT) THEN
              IF    (KDR .GE. IDRT) THEN
C--------------- i.e., IDRT=1 or 2
                XFINAL = XLIM
                GOTO 6
              ELSEIF(KDR .LT. IDRT) THEN
C--------------- i.e., IDRT > 2
C---------------- S'il y a (IDRT-1)>2 droites en sorties (e.g., plans 
C                                    de detecteurs)
                KDR = KDR + 1
                CT=COS(T)
                ST=SIN(T)
                GOTO 1
              ENDIF
            ELSEIF(MIRROR) THEN
              IF(BACK) THEN
                PAF = -PAF
                CALL DEPLA(PAF)
                CALL COFIN
     >          (KART,NPLT,LST,PAF,Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT,EVNT,
     >             *97)
C     >             *99)
                IF(LST .EQ. 2) THEN 
                  CALL IMPPLA(NPLT,PAF,Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT)
                ELSEIF(LST .EQ. 3) THEN 
                  CALL IMPPLA(  30,PAF,Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT)
                ENDIF
                XFINAL = -CM(IDRT)
                GOTO 6
              ENDIF
            ENDIF
C------------- if .not.mirror etc. 

          ENDIF
        ENDIF
      ENDIF

C-----------------------    TEST  EFB's in Multipoles etc.
C                              IFB set to .ne.0 due to mixed sharp edge+FF
      IF(IFB .NE. 0 ) THEN
C------- Entrance or exit EFB's 

        IF( KART .EQ. 1) THEN
C--------- COORDONNEES CARTESIENNES

          KFB = 0

          IF( AFB(1)*X + BFB(1)*Y + CFB(1) .LE. 0.D0 ) THEN

            IF( IFB .EQ. -1 .OR. IFB .EQ. 2 ) KFB=1
C----------- EFB entree

          ELSEIF(AFB(2)*X + BFB(2)*Y + CFB(2) .LE. 0.D0 ) THEN

            IF( IFB .EQ. 1 .OR. IFB .EQ. 2 ) KFB=2
C----------- EFB sortie
          ENDIF

          IF( KFB .NE. 0 ) THEN
            D = ABS(AFB(KFB)*X + BFB(KFB)*Y + CFB(KFB))
            IF(D .LT. DXI) THEN
              AL = AFB(KFB)
              BL = BFB(KFB)
              CL = -D
              ST = SIN(T)
              CT = COS(T)
              COSTA = AL*CT + BL*ST
              PAF = D/COSTA
              CALL DEPLA(PAF)
              COSTA=AL*CT+BL*ST
              CALL ITER(AL,BL,CL,PAF,COSTA,KEX,*97)
C              CALL ITER(AL,BL,CL,PAF,COSTA,KEX,*99)
              CALL COFIN
     >          (KART,NPLT,LST,PAF,Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT,EVNT,
     >           *97)
C     >           *99)
              CT=COS(T)
              ST=SIN(T)
              IF(ABS(D) .GT. DLIM) GOTO 1
            ENDIF
          ENDIF
C------------ kfb .ne. 0

        ENDIF
C---------- kart .eq. 1

      ENDIF
C-------- ifb .ne. 0

      DX=XLIM-X
C        write(89,*)  x, dxi, dx, xlim, ' 1  sbr integr' 
      IF(DX/DXI .LT. 1.D0) THEN
        IF(.NOT.MIRROR) THEN
          GOTO 2
        ELSE
          IF(.NOT.BACK) CALL KSTOP(8,IT,KEX,*97)
        ENDIF
      ENDIF
      PAF = PAS
      CALL DEPLA(PAF)
      CALL COFIN
     >  (KART,NPLT,LST,PAF,Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT,EVNT,
     >   *97)
C     >   *99)

C rustine RACCAM -----------
c        if(it.eq.1.and.noel.eq.8) then
c          if(x.gt.-.0046) then 
c            if(first) then
c              write(*,*) ' *****it,x,t,noel,first : ', it,x,t,noel,first
c              t = T - 15.7e-3
c              first = .false. 
c              write(*,*) ' *****it,x,t,noel,first : ', it,x,t,noel,first
c            endif
c          endif
c        endif
C---------------------------

      CT=COS(T)
      ST=SIN(T)

      X2 = X
      Y2 = Y
      Z2 = Z

      GOTO 1
C---------  FIN DE BOUCLE SUR DXI
 
 2    CONTINUE
C        write(89,*)  x, dxi, dx, xlim, ' 2  sbr integr' 
      IF(KART .EQ. 1) THEN
        AL=1.D0
        BL=0.D0
        CL=-DX
        PAF=DX/CT
      ELSEIF(KART .EQ. 2) THEN
        AL=COS(DX)
        BL=-SIN(DX)
        CL=BL*Y
        PAF=-CL/CT
      ENDIF
 
      CALL DEPLA(PAF)
      COSTA=AL*CT+BL*ST
      CALL ITER(AL,BL,CL,PAF,COSTA,KEX,*97)
      CALL COFIN
     >  (KART,NPLT,LST,PAF,Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT,EVNT,
     >   *97)
C        write(89,*)  x, dxi, dx, xlim, ' 3  sbr integr' 

      XFINAL = XLIM

 6    CONTINUE
      DX=XFINAL-X

      IF(ABS(DX) .GT. DLIM) THEN
        CT=COS(T)
        ST=SIN(T)
        X=XFINAL
        Y=Y+(DX*ST)/CT
        Z=Z+(DX*TAN(P))*(1.D0/CT)
        PAF = DX/(CT*COS(P))
        SAR= SAR+PAF
        IF(QT*AMT .NE. 0.D0) THEN 
          QBRO = BR*CL9*QT
          DTAR = PAF / (QBRO/SQRT(QBRO*QBRO+AMT*AMT)*CL9)
          TAR = TAR + DTAR
        ENDIF
        IF(EVNT) THEN
          IF(KSPN .EQ. 1) THEN
C----------- Spin tracking
            IF(KART .EQ. 2) CALL SPNROT(IT,ZERO,ZERO,-DX)
            CALL SPNTRK(IT,PAF)
          ELSE         
            CALL EVENT(
     >      PAF,Y,T,Z,P,X,-DX,ZERO,ZERO,ZERO,BR,SAR,TAR,KEX,IT,
     >      AMT,QT,BORO,KART,KSPN,IFDES,KGA,KSYN,IMAX,*97)
          ENDIF
        ENDIF
      ENDIF
C        write(89,*)  x, dxi, dx, xlim, ' 4  sbr integr' 

 99   CONTINUE
      IF(KALC .EQ. 2 ) THEN
        CALL CHAMK(X,Y,Z,*97)
      ELSE

CCCCCCCCCCCCC FM, rustine Synchrotron Radiation LHC / Zpop  CCCC
C *.99999D0 introduced to avoid that B be set to zero 
C in SBR BEND due to the  test X.LT.XS in sharp edge magnet case
        CALL CHAMC(X*.99999D0,Y,Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        CALL CHAMC(X,Y,Z)
CCCCCCCCCCCCCC

      ENDIF

C-------- Wedge correction in BEND, in MULTIPOL with non zero B1, etc.
      IF(WEDGS) CALL WEDGKI(2,T,Z,P,WDGS,FFXTS)

 97   CONTINUE

C----- Print last step if requested
      IF(LST .GE. 1) THEN
        IF(LST .EQ. 2) THEN
          CALL IMPPLA(NPLT,PAF,Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT)
        ELSE
          CALL IMPDEV
        ENDIF
      ENDIF

      RETURN

      ENTRY INTEG1(WDGI,FFXTEI)
      WEDGE = .TRUE.
      WDGE=WDGI
      FFXTE = FFXTEI
      RETURN
      ENTRY INTEG2(WDGI,FFXTSI)
      WEDGS = .TRUE.
      WDGS = WDGI
      FFXTS = FFXTSI
      RETURN
      ENTRY INTEG3
      WEDGE = .FALSE.
      WEDGS = .FALSE.
      RETURN
      ENTRY INTEG4(NSTP)
      NSTP = NSTEP
      RETURN
      ENTRY INTEG5(
     >             PASO)
      PASO=PAS
      RETURN
      ENTRY INTEG6(AREGI,BREGI,CREGI)
      AREG(1)=AREGI(1)
      BREG(1)=BREGI(1)
      CREG(1)=CREGI(1)
      AREG(2)=AREGI(2)
      BREG(2)=BREGI(2)
      CREG(2)=CREGI(2)
      RETURN
      END
