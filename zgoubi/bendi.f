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
      SUBROUTINE BENDI(SCAL,
     >                      XL,DEV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.AIM.H"     ! COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      INCLUDE "C.CHAFUI.H"     ! COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.DROITE.H"     ! COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      PARAMETER(MPOL=10)
      INCLUDE "C.MULTPL.H"   ! COMMON/MULTPL/ BM(MPOL),DLE(MPOL),DLS(MPOL),DI(MPOL,MCOEF),DS(MPOL,MCOEF),RTB(MPOL)
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE "C.SYNRA.H"     ! COMMON/SYNRA/ KSYN

      EQUIVALENCE (RTB(1),CTE),(RTB(2),STE),(RTB(4),CTS),(RTB(5),STS)

      LOGICAL SHARPE, SHARPS
      DIMENSION  AREG(2),BREG(2),CREG(2)
      SAVE DEVO

C----- Magnet length, geometrical
      XL =A(NOEL,10)
C----- Skew angle
      BM(6) = A(NOEL,11)
C----- Field
      BM(1)  =A(NOEL,12)*SCAL
      GAP  =A(NOEL,13)
      AK1    =A(NOEL,14)*SCAL
      AK2   =A(NOEL,15)*SCAL
      IF(GAP .NE. 0.D0) THEN
        G1 = AK1/GAP
        G2 = AK2/GAP**2
      ELSE
        G1 = 0.D0
        G2 = 0.D0
      ENDIF

      DEV = 2.D0 * ASIN(XL/2.D0/(BORO/BM(1)))
      KP = NINT(A(NOEL,70))

c      write(88,*) ' bendi A(NOEL,73) ', A(NOEL,73)


      IF( KP .EQ. 3 ) THEN
         IF(A(NOEL,73) .NE. 0.D0) DEV = -A(NOEL,73) * 2.D0
      ENDIF

c         write(88,*) '  ',2.D0 * ASIN(XL/2.D0/(BORO/BM(1))),
c     > A(NOEL,73), -A(NOEL,73) * 2.D0,XL,BORO,BM(1),scal

      IF(BM(6) .NE. 0.D0) THEN
C------- BEND is skewed
        DEVH=ATAN(TAN(DEV)*COS(BM(6)))
        DEVV=ATAN(TAN(DEV)*SIN(BM(6)))
      ENDIF

      XE = A(NOEL,20)
      DLE(1) = A(NOEL,21)
      WE = A(NOEL,22)
      TE=.5D0*DEV-WE
      XLS = A(NOEL,40)
      DLS(1) = A(NOEL,41)
      WS = A(NOEL,42)
      TS=.5D0*DEV-WS
      DO 8 I=1,6
        CE(I) =A(NOEL,30+I)
        CS(I) =A(NOEL,50+I)
 8    CONTINUE
      IOPIN = NINT(A(NOEL,30))
      DISTE = -A(NOEL,30)
      IOPOU = NINT(A(NOEL,50))
      DISTS = -A(NOEL,50)
      CALL  BENDFI(IOPIN,DISTE,IOPOU,DISTS)

C----------- Champ DE FUITE
      SHARPE=DLE(1) .LE. 0.D0
      SHARPS=DLS(1) .LE. 0.D0

      IF(SHARPE) THEN
        FINTE = XE
C        XE=0.D0
C        Introduced so to be consistent with xls below, yet actually not necessary...
C        if (xe.lt.1d-10) xe = 5.d0*tan(abs(te))
        GAPE = -DLE(1)
      ENDIF
      IF(SHARPS) THEN
        FINTS = XLS
        XLS=0.D0
C        In INTEGR "Droite" has to be intersected by trajectory within X<XLIM (X cannot
C        be > XLIM), hence introducing some additional integration distance beyond XLIM
C        so to encompass "Droite".
C FM, 130605. Removed, causes problems (see zgoubi/folks/samTygier/problem...)
C FM 171122. Coded step size with wedge BEND causes problems,
C not compatible with "3 regions" (it remains to determine why !),
C use cm step instead  (see [pathTo]/SVN/current/debugging/spinFlipper)
        IF (XLS.LT.1D-10) THEN
          CALL CHXC1R(
     >                KPAS)
          IF(KPAS .EQ. 0) THEN
            XLS = 5.D0*TAN(ABS(TS))
          ELSE
            CALL ENDJOB(' Pgm bendi. Entrance/body/exit step size mode '
     >      //'is not compatible with EXIT WEDGE ANGLE.'
     >      //' Instead, use explicit single step size value. ',-99)
          ENDIF
        ENDIF
        GAPS = -DLS(1)
      ENDIF
      XI = 0.D0
      XLIM = XL + XE + XLS
      XF = XLIM
      XS = XL + XE
      CTE=COS(TE)
      STE=SIN(TE)
      CTS=COS(TS)
      STS=SIN(TS)

C----- SHARP EDGE => INTEGR STOPPE SUR DR. DE COUPURE
      IF(SHARPE) THEN
C------- Correction for entrance wedge
        CALL INTEG1(-TE,FINTE,GAPE)

        IDRT = -1
        CA(1)=CTE
        SA(1)=STE
        CM(1)=-XE*CA(1)
      ELSE
        DI(1,1)=-BM(1)/DLE(1)
        DO 10 I=2,4
          DI(1,I)=-DI(1,I-1)/DLE(1)
 10     CONTINUE
      ENDIF

      IF(SHARPS) THEN
C------- Correction for exit wedge
        CALL INTEG2(-TS,FINTS,GAPS)

        IF(IDRT .EQ. -1) THEN
          IDRT = 2
        ELSE
          IDRT = 1
        ENDIF
        CA(2)= CTS
        SA(2)=-STS
        CM(2)=-XS*CA(2)
      ELSE
        DS(1,1)=-BM(1)/DLS(1)
        DO 11 I=2,4
          DS(1,I)=-DS(1,I-1)/DLS(1)
 11     CONTINUE
      ENDIF

      CALL CHXC1R(
     >            KPAS)
      IF(KPAS.GE.1) THEN
        IF(XE+XLS .GE.XL) THEN
          CALL ENDJOB(' Pgm bendi. Entrance/body/exit step size mode '
     >    //'is not compatible with fringe fields overlapping in body.'
     >    //' Instead, use explicit single step size value. ',-99)
        ENDIF
        AREG(1)=CTE
        BREG(1)=STE
        CREG(1)=-2.D0*XE*AREG(1)
        AREG(2)=CTS
        BREG(2)=-STS
        CREG(2)=(-2.D0*XS+XLIM)*AREG(2)
        CALL INTEG6(AREG,BREG,CREG)
      ENDIF

      IF(NRES.GT.0) THEN
        WRITE(NRES,100) ' BEND',XL,BORO/BM(1)*DEV,DEV*DEG,DEV
     >       , GAP, G1, G2
 100    FORMAT(1P, /,5X,' +++++  ',A10,'  : ',
     >       //,15X, ' Length    = ',E14.6,' cm'
     >       ,/,15X, ' Arc length    = ',E14.6,' cm'
     >       ,/,15X, ' Deviation    = ',E14.6,' deg.,  ',E14.6,' rad'
     >       ,/,15X, ' GAP   = ',E14.6,' cm'
     >       ,/,15X, ' Gradient   = ',E14.6,' kG/cm'
     >       ,/,15X, ' Grad-prime   = ',E14.6,' kG/cm^2',/)
        WRITE(NRES,103) BM(1),BM(1)/SCAL,BORO/BM(1)
 103    FORMAT(1P,15X,' Field  =',E15.7,'  kG ',
     >  '  (i.e., ',E15.7,' * SCAL)',
     >  /,15X, ' Reference curvature radius (Brho/B) = ',E15.7,' cm')
        WRITE(NRES,105) BM(6)
 105    FORMAT(1P,15X, ' Skew  angle  = ',E14.6,'  rad')
        IF(BM(6).NE.0.D0) WRITE(NRES,FMT='(15X,
     >       '' Projected  deviations  in H/V  planes  = '',
     >               1P,G13.6,''/'',G13.6,''  (rad)'')') DEVH, DEVV

C        WRITE(NRES,104) 'D''ENTREE'
C 104    FORMAT(/,15X,' FACE  ',A)
        WRITE(NRES,104) 'Entrance '
 104    FORMAT(/,15X,A9,' face  ')
        WRITE(NRES,101) XE,DLE(1),WE
 101    FORMAT(15X,' DX = ',F10.3,'    LAMBDA = ',F10.3
     >        ,/,15X,' Wedge  angle  =',F10.6,' RD')
C     >        ,/,15X,' ANGLE  DE  COIN =',F10.6,' RD')
        IF( .NOT. SHARPE ) WRITE(NRES,132) (CE(I),I=1,6)
C 132    FORMAT(15X,' COEFFICIENTS DE Champ DE FUITE:',/,16X,6F9.5)
 132    FORMAT(15X,' Fringe  field  coefficients :',/,16X,6F9.5)

C        WRITE(NRES,104) 'DE  SORTIE'
        WRITE(NRES,104) 'Exit     '
        WRITE(NRES,101) XLS,DLS(1),WS
        IF( .NOT. SHARPS ) WRITE(NRES,132) (CS(I),I=1,6)

        IF( SHARPE)
     >    WRITE(NRES,FMT='(/,''  ***  Warning : entrance sharp edge '',
     >    ''entails vertical wedge focusing approximated with'',
     >    '' first order kick, FINT values entr/exit : '',1P,G12.4)')
     >         FINTE
        IF( SHARPS)
     >    WRITE(NRES,FMT='(/,''  ***  Warning : exit sharp edge '',
     >    ''entails vertical wedge focusing approximated with'',
     >    '' first order kick, FINT values entr/exit : '',1P,G12.4)')
     >        FINTS
      ENDIF

      CALL BENDF2(G1,G2)

      DEVO = DEV

C----- SR-LOSS SWITCHED ON BY PROCEDURE SRLOSS
      IF(KSYN.GE.1) THEN
        IF(KFLD .EQ. MG) THEN
          IF(BM(1).NE.0.D0) CALL SYNPAR(BM(1),BM(6),XL)
        ENDIF
      ENDIF

      RETURN

      ENTRY BENKL(
     >            DEVOO)
      DEVOO = DEVO
      RETURN
      END
