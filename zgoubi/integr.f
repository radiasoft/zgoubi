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
      SUBROUTINE INTEGR(EVNT,BACKW,MIRROR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EVNT,BACKW,MIRROR
C     -------------------------------------------------------------
C     APPELE PAR TRANSF.
C     CALCULE UNE TRAJECTOIRE DE X A XLIM; DXI=PAS EN X OU EN ANGLE
C     PAS = PAS EN S : cf. DEPLA, COFIN
C     -------------------------------------------------------------
      INCLUDE "C.AIM.H"     ! COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      INCLUDE "C.CHAFUI.H"     ! COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
      INCLUDE "C.CHAVE_2.H"     ! COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      INCLUDE "C.CONST_3.H"      ! COMMON/CONST/ CL9,CEL,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "C.DROITE.H"     ! COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      INCLUDE "C.EFBS.H"     ! COMMON/EFBS/ AFB(2), BFB(2), CFB(2), IFB
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE "C.MARK.H"     ! COMMON/MARK/ KART,KALC,KERK,KUASEX
      INCLUDE "MAXTRA.H"
C      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
CC     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
C      INCLUDE "C.GASC.H"     ! COMMON/GASC/ AI, DEN, KGA
      INCLUDE "MAXCOO.H"
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INCLUDE "C.ORDRES.H"     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
C FM Dec.2010, KSPN needed at PAF
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      INCLUDE "C.STEP.H"     ! COMMON/STEP/ TPAS(3), KPAS
C      INCLUDE "C.SYNRA.H"     ! COMMON/SYNRA/ KSYN
      INCLUDE "C.TRAJ.H"     ! COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT

Ccccccccccccc  SPACE CHARGE MALEK cccccccccccc
C      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT)
      logical TSPCH
C      COMMON/SPACECHA/ TSPCH,TLAMBDA,Rbeam,Xave,Emitt,Tave,Bunch_len,
      COMMON/SPACECHA/ TLAMBDA,Rbeam(2),Xave(2),Emitt(2),Tave,Bunch_len,
     >                Emittz, BTAG, SCkx, SCky, TSPCH
cc      INCLUDE 'MXLD.H'              ! pour faire les tests: to be deleted please 
cc      COMMON /LABEL/ LABEL(MXL,2)   ! pour faire les tests: to be deleted please
C       INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
cccccccccccccccccccccccccccccccccccccccccccccc

      LOGICAL WEDGE, WEDGS
      SAVE WEDGE, WEDGS
 
      PARAMETER (DLIM = 1.D-6)
      INCLUDE 'MXSTEP.H'
      PARAMETER (PI24= 2.467401101D0)
 
      LOGICAL BACK, CHREG
 
      SAVE NSTEP
      DIMENSION  AREG(2),BREG(2),CREG(2), AREGI(2),BREGI(2),CREGI(2)
      SAVE AREG,BREG,CREG,KREG
 
      SAVE WDGE, WDGS, FINTE, FINTS, GAPE, GAPS
 
      LOGICAL FITMEM
      SAVE FITMEM
 
      LOGICAL CONSTY, CSTYI, IDLUNI
      SAVE CONSTY 

      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KEY

      DIMENSION F(MXJ,MXT)
 
      DATA WEDGE, WEDGS / .FALSE.,  .FALSE./
      DATA WDGE, WDGS, FINTE, FINTS, GAPE, GAPS / 6*0.D0 /
 
      DATA FITMEM / .FALSE. /
      DATA CONSTY / .FALSE. /

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
      IF(FITMEM) CALL FITMM2(IT)
 
C----- DEBUT DE BOUCLE SUR DXI
C      Start loop on DXI
    1 CONTINUE
 
C---------------------  Some tests to possibly stop integration
      NSTEP = NSTEP+1
      IF(NSTEP .GT. MXSTEP) THEN
        CALL ZGNOEL(
     >             NOEL)
        CALL ZGKLEY(  
     >              KEY)
        WRITE(ABS(NRES),fmt='(a,2(i0,a))') 
     >  'Particle exceeded maximum # of steps allowed : ',MXSTEP,
     >  '.  At element number ',NOEL,' ('//KEY//')'
        WRITE(6,fmt='(a,2(i0,a))') 
     >  'Particle exceeded maximum # of steps allowed : ',MXSTEP,
     >  '.  At element number ',NOEL,' ('//KEY//')'
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
          IF(NRES.GT.0) WRITE(NRES,100)Y,X,T,IT,' went out of'
  100     FORMAT (/,4X,'Y =',F8.2,' X =',F8.2,' T =',F6.3
     >    ,' Trajectory # ',I7,A,' field region.')
          ISORT=2
          GOTO 7
        ENDIF
        IF(ISORT .EQ. 1)  GOTO 7
        IF(NRES.GT.0) WRITE(NRES,100) Y,X,T,IT,' came back in.'
        ISORT=1
 
      ELSE
C-------  Compute  B(X,Y,Z), E(X,Y,Z)  from  mathematical  2D or 3D  field  models
 
        CALL CHAMC(X2,Y2,Z2)
 
      ENDIF
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C------- Correction for space charge calculation --Begin--
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          IF (TSPCH)  THEN


ccccccccccccc   get the path legnth from the first central orbit   cccccccccccccc
           if (IT .EQ. 1) then
                   tlength0=tlength0+PAS
                   tlength=tlength0 !tlength is the element path length in cm
                   if (tlength .LE. 1.1*PAS) then
                    T1=T
                    Y1=Y
                    Z1=Z
                   endif 
           elseif (IT .NE. 1) then                   
                   tlength0=0.0
           endif        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                


           if ((IT .NE. 1) .and. (IT .NE. ITinit)) then
             Yinit=Y
             Zinit=Z
             ITinit=IT
! correct             Tsp=SCkx*(Yinit-Y1)*tlength
! correct             Psp=SCky*(Zinit-Z1)*tlength
             Tsp=SCkx*(Yinit-Xave(1)*100.0)*tlength
             Psp=SCky*(Zinit-Xave(2)*100.0)*tlength

             T=T+Tsp   
             P=P+Psp

           endif   

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccc    write the data to a test4 file   ccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c           if ((IT .EQ. 1) .and. (tlength .LE. 2*PAS)) then
c             call ZGIPAS( 
c     >                   IPASS,NRBLT)
c              write(88,fmt=
c     >        '(1p,9(e14.6,1x),i10,1x,e14.6,a)') 
c     >        Rbeam(1),Rbeam(2),Xave(1),SAR,Emitt(1),Emitt(2),Tave,
c     >        Bunch_len, Emittz,IPASS,BTAG,' integr spach'
c           endif  
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  
          ENDIF   


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C------- Correction for space charge calculation --End--
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


 7    CONTINUE
 
      IF(NSTEP .EQ. 1) THEN
C-------- Entrance wedge correction in BEND, in MULTIPOL(if non zero B1),
        IF(WEDGE) CALL WEDGKI(1,T,Z,P,WDGE,FINTE,GAPE)
      ENDIF

      CALL DEVTRA
 
      IF(KPAS.NE.0 ) THEN
C------- Test case diffrnt step values @ entrance|body|exit.
C         Exist entrance and/or exit fringe field regions and XPAS is coded.
 
C-------- CHREG is .true. if particle is going to next region
        IF(CHREG(KART,X,Y,AREG,BREG,CREG,DXI,TPAS,
     >                                            KREG,AL,BL,D)) THEN
            CL = -D
            ST = SIN(T)
            CT = COS(T)
            COSTA = AL*CT + BL*ST
            PAF = D/COSTA
            CALL DEPLA(PAF)
            COSTA=AL*CT+BL*ST
            CALL ITER(AL,BL,CL,PAF,COSTA,KEX,IT,*97)
            CALL COFIN(KART,NPLT,LST,PAF,KEX,IT,AMT,QT,EVNT,
     >                                            Y,T,Z,P,X,SAR,TAR,*97)
            CT=COS(T)
            ST=SIN(T)
            PAS = TPAS(KREG)
            DXI = PAS
 
            GOTO 1
        ENDIF
      ENDIF
 
C FM 08/99      IF(IDRT .GE. 2) THEN
      IF(IDRT .GE. 1) THEN
C------- DROITE(S) DE COUPURE EN SORTIE
C        "Droite" has to be intersected by trajectory within X<XLIM (X cannot
C        be > XLIM)
 
        IF    (KART .EQ. 1) THEN
C--------- COORDONNEES CARTESIENNES
          D = ABS(CA(KDR)*X + SA(KDR)*Y + CM(KDR))
 
          IF(D .LT. DXI) THEN
            AL = CA(KDR)
            BL = SA(KDR)
            CL = -D
            ST = SIN(T)
            CT = COS(T)
            COSTA = AL*CT + BL*ST
            PAF = D/COSTA
 
            IF(.NOT.MIRROR) THEN
C------------- Step onto the exit 'DROITE'
 
              CALL DEPLA(PAF)
              COSTA=AL*CT+BL*ST
              CALL ITER(AL,BL,CL,PAF,COSTA,KEX,IT,*97)
              CALL COFIN(KART,NPLT,LST,PAF,KEX,IT,AMT,QT,EVNT,
     >                                            Y,T,Z,P,X,SAR,TAR,*97)
              IF(LST .EQ. 2) THEN
                CALL IMPPLA(NPLT,Y,T,Z,P,X,SAR,TAR,PAF,AMT,QT,KEX,IT)
              ELSEIF(LST .EQ. 3) THEN
                CALL IMPPLA(  30,Y,T,Z,P,X,SAR,TAR,PAF,AMT,QT,KEX,IT)
              ENDIF
 
C FM 08/99              IF    (KDR .EQ. IDRT) THEN
              IF    (KDR .GE. IDRT) THEN
C--------------- i.e., IDRT=1 or 2
                XFINAL = XLIM
                GOTO 6
              ELSEIF(KDR .LT. IDRT) THEN
C---------------- Thus, IDRT > 2
C---------------- There are (IDRT-1)>2 "droites en sorties" (e.g., spectrometer detector planes)
                KDR = KDR + 1
                CT=COS(T)
                ST=SIN(T)
                GOTO 1
              ENDIF
            ELSEIF(MIRROR) THEN
              IF(BACK) THEN
                PAF = -PAF
                CALL DEPLA(PAF)
                CALL COFIN(KART,NPLT,LST,PAF,KEX,IT,AMT,QT,EVNT,
     >                                            Y,T,Z,P,X,SAR,TAR,*97)
                IF(LST .EQ. 2) THEN
                  CALL IMPPLA(NPLT,Y,T,Z,P,X,SAR,TAR,PAF,AMT,QT,KEX,IT)
                  CALL IMPPLA(  30,Y,T,Z,P,X,SAR,TAR,PAF,AMT,QT,KEX,IT)
                ELSEIF(LST .EQ. 3) THEN
                  CALL IMPPLA(  30,Y,T,Z,P,X,SAR,TAR,PAF,AMT,QT,KEX,IT)
                ENDIF
                XFINAL = -CM(IDRT)
                GOTO 6
              ENDIF
            ENDIF
C------------- if .not.mirror etc.
 
          ENDIF
 
        ELSEIF( KART .EQ. 2) THEN
C--------- Cylindrical coordinates
C Array CA contains angle of integration boundary
 
          CALL ENDJOB(
     >    ' Sbr integr : IDRT not implemented for polar frame',-99)
 
          IF(ABS(CA(2)-X) .LT. DXI) THEN
 
          ENDIF ! XLIM-X
 
        ENDIF  ! KART
      ENDIF ! IDRT .GE. 1
 
 
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
              CALL ITER(AL,BL,CL,PAF,COSTA,KEX,IT,*97)
              CALL COFIN(KART,NPLT,LST,PAF,KEX,IT,AMT,QT,EVNT,
     >                                            Y,T,Z,P,X,SAR,TAR,*97)
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

      IF(DX/DXI .LT. 1.D0) THEN
        IF(.NOT.MIRROR) THEN
          GOTO 2
        ELSE
          IF(.NOT.BACK) CALL KSTOP(8,IT,KEX,*97)
        ENDIF
      ENDIF
 
      PAF = PAS
      CALL DEPLA(PAF)
 
      CALL COFIN(KART,NPLT,LST,PAF,KEX,IT,AMT,QT,EVNT,
     >                                            Y,T,Z,P,X,SAR,TAR,*97)
 
      IF(FITMEM) CALL FITMM(IT,Y,T,Z,P,SAR,DP,TAR,PAS)
 
      CT=COS(T)
      ST=SIN(T)
 
C  A trick for tests at constant coordinate -----------------
         IF(CONSTY) THEN
               CALL MAJTR1(
     >                     F)
                        Y = F(2,IT)
                        Z = F(4,IT)
                     T = 0.D0
                     P = 0.D0  
                    X = NSTEP*DXI
         ENDIF
C------------------------------------------------------------
 
      X2 = X
      Y2 = Y
      Z2 = Z
 
      GOTO 1
C---------  FIN DE BOUCLE SUR DXI
C           End loop on DXI
 
 2    CONTINUE
 
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
      CALL ITER(AL,BL,CL,PAF,COSTA,KEX,IT,*97)
      CALL COFIN(KART,NPLT,LST,PAF,KEX,IT,AMT,QT,EVNT,
     >                                            Y,T,Z,P,X,SAR,TAR,*97)
 
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
        QBRO = QBR*CL9
        DTAR = PAF / (QBRO/SQRT(QBRO*QBRO+AMT*AMT)*CL9)
        TAR = TAR + DTAR
        IF(EVNT) THEN
          IF(KSPN .EQ. 1) THEN
            IF(IDRT .EQ. 0) THEN   ! FM. Added Jan 2017. Pb was : spntrk won't linearly extrapolate as B .ne. 0
C----------- Spin tracking
              IF(KART .EQ. 2) CALL SPNROT(IT,ZERO,ZERO,-DX)
C              CALL SPNTRK(IT,PAF)
              CALL SPNTRK(PAF)
            ENDIF
          ENDIF
c          ELSE

          CALL EVENT(
     >    PAF,Y,T,Z,P,X,ZERO,QBR,SAR,TAR,KEX,IT,
     >    AMT,Q,BORO,KART,IMAX,*97)

        ENDIF
      ENDIF
 
      IF(FITMEM) CALL FITMM(IT,Y,T,Z,P,SAR,DP,TAR,PAS)
 
 99   CONTINUE
c      IF(KALC .EQ. 2 ) THEN
c        CALL CHAMK(X,Y,Z,*97)
c      ELSE
CCCCCCCCCCCCC FM, rustine Synchrotron Radiation LHC / Zpop  CCCC
C *.999999999D0 introduced to avoid that B be set to zero
C in SBR BEND due to the  test X.LT.XS in sharp edge magnet case
C        CALL CHAMC(X*.999999999D0,Y,Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c        CALL CHAMC(X,Y,Z)
CCCCCCCCCCCCCC
c      ENDIF
 
C-------- Wedge correction in BEND, in MULTIPOL with non zero B1, etc.
      IF(WEDGS) CALL WEDGKI(2,T,Z,P,WDGS,FINTS,GAPS)
 
 97   CONTINUE
 
C----- Print last step
      IF(LST .GE. 1) THEN
        IF(LST .EQ. 2) THEN
          CALL IMPPLA(NPLT,Y,T,Z,P,X,SAR,TAR,PAF,AMT,QT,KEX,IT)
        ELSE
          CALL IMPDEV
        ENDIF
      ENDIF
 
      RETURN
 
      ENTRY INTEG1(WDGI,FINTEI,GAPEI)
      WEDGE = .TRUE.
      WDGE=WDGI
      FINTE = FINTEI
      GAPE = GAPEI
      RETURN
 
      ENTRY INTEG2(WDGI,FINTSI,GAPSI)
      WEDGS = .TRUE.
      WDGS = WDGI
      FINTS = FINTSI
      GAPS = GAPSI
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
 
      ENTRY INTEG8(KFITMM)
      FITMEM = KFITMM.EQ.1
      RETURN

      ENTRY INTEGA(CSTYI)
      CONSTY = CSTYI
      RETURN
 
      END
