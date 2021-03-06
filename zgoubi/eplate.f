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
      SUBROUTINE EPLATE(MPOL,SCAL,
     >          DEV,XL,EM,DLE,DLS,DE,DS,XE,XS,CE,CS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EM(*),DLE(*),DLS(*),DE(MPOL,*),DS(MPOL,*)
      PARAMETER(MCOEF=6)
      DIMENSION CE(MCOEF), CS(MCOEF)

      INCLUDE "C.AIM.H"     ! COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.DROITE.H"     ! COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE "C.SYNRA.H"     ! COMMON/SYNRA/ KSYN

C      EQUIVALENCE (RTB(1),CTE),(RTB(2),STE),(RTB(4),CTS),(RTB(5),STS)

      LOGICAL SHARPE, SHARPS
      DIMENSION  AREG(2),BREG(2),CREG(2)

C----- Length
      XL =A(NOEL,10)
C----- Skew angle
      EM(6) = A(NOEL,11)
C----- V1, V2, d
      EM(1)  =A(NOEL,12)*SCAL
      EM(2)  =A(NOEL,12)*SCAL
      gap  =A(NOEL,12)*SCAL

      DEV = 2.D0 * ASIN(XL/2.D0/((BORO*(DPREF+HDPRF))/EM(1)))
      KP = NINT(A(NOEL,70))
      IF( KP .EQ. 3 ) THEN
        IF(A(NOEL,73) .NE. 0.D0) DEV = -A(NOEL,73) * 2.D0
      ENDIF

      IF(EM(6) .NE. 0.D0) THEN
C------- BEND is skewed
        DEVH=ATAN(TAN(DEV)*COS(EM(6)))
        DEVV=ATAN(TAN(DEV)*SIN(EM(6)))
      ENDIF

      XE = A(NOEL,20)
      DLE(1) = A(NOEL,21)
      WE = A(NOEL,22)
      TE=0.D0  !.5D0*DEV-WE
      XLS = A(NOEL,40)
      DLS(1) = A(NOEL,41)
      WS = A(NOEL,42)
      TS=0.D0 !!.5D0*DEV-WS
      DO 8 I=1,6
        CE(I) =A(NOEL,30+I)
        CS(I) =A(NOEL,50+I)
 8    CONTINUE

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
        if (xls.lt.1d-10) xls = 5.d0*tan(abs(ts))
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
        IDRT = -1
        CA(1)=CTE
        SA(1)=STE
        CM(1)=-XE*CA(1)
      ELSE
        DE(1,1)=-EM(1)/DLE(1)
        DO 10 I=2,4
          DE(1,I)=-DE(1,I-1)/DLE(1)
 10     CONTINUE
      ENDIF

      IF(SHARPS) THEN
        IF(IDRT .EQ. -1) THEN
          IDRT = 2
        ELSE
          IDRT = 1
        ENDIF
        CA(2)= CTS
        SA(2)=-STS
        CM(2)=-XS*CA(2)
      ELSE
        DS(1,1)=-EM(1)/DLS(1)
        DO 11 I=2,4
          DS(1,I)=-DS(1,I-1)/DLS(1)
 11     CONTINUE
      ENDIF

      CALL CHXC1R(
     >            KPAS)
      IF(KPAS.GE.1) THEN
        IF(XE+XLS .GE.XL) THEN
          CALL ENDJOB(' Pgm eplate. Entrance/body/exit step size mode '
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
        WRITE(NRES,100) ' EPLATES',XL,EM(1),EM(2),GAP
 100    FORMAT(1P, /,5X,' +++++  ',A10,'  : ',
     >       //,15X, ' Length    = ',E14.6,' cm'
     >       ,/,15X, ' V1    = ',E14.6,' V'
     >       ,/,15X, ' V1    = ',E14.6,' V'
     >       ,/,15X, ' Gap    = ',E14.6,' cm ',/)
        WRITE(NRES,105) EM(6)
 105    FORMAT(1P,15X, ' Skew  angle  = ',E14.6,'  rad')
        IF(EM(6).NE.0.D0) WRITE(NRES,FMT='(15X,
     >       '' Projected  deviations  in H/V  planes  = '',
     >               1P,G13.6,''/'',G13.6,''  (rad)'')') DEVH, DEVV

        WRITE(NRES,104) 'Entrance '
 104    FORMAT(/,15X,A9,' face  ')
        WRITE(NRES,101) XE,DLE(1),WE
 101    FORMAT(15X,' DX = ',F10.3,'    LAMBDA = ',F10.3
     >        ,/,15X,' Wedge  angle  =',F10.6,' RD')
        IF( .NOT. SHARPE ) WRITE(NRES,132) (CE(I),I=1,6)
 132    FORMAT(15X,' Fringe  field  coefficients :',/,16X,6F9.5)

        WRITE(NRES,104) 'Exit     '
        WRITE(NRES,101) XLS,DLS(1),WS
        IF( .NOT. SHARPS ) WRITE(NRES,132) (CS(I),I=1,6)

        IF( SHARPE .OR. SHARPS )
     >    WRITE(NRES,FMT='(/,''  ***  Warning : sharp edge '',
     >    ''entails non-conservation of energy.'')')

      ENDIF

      RETURN
      END
