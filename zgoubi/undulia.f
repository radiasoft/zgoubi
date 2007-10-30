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
      SUBROUTINE UNDULIA(SCAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
C      COMMON/CHAFUI/ XE,XS,CE(6),CS(6),QCE(6),QCS(6)
      COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      PARAMETER(MPOL=10)
      COMMON/MULTPL/ BBM(MPOL),DLE(MPOL),DLS(MPOL)
     >,DFM(MPOL),DE(MPOL,10),DS(MPOL,10),RTB(MPOL)
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,BR
 
      EQUIVALENCE (RTB(1),CTE),(RTB(2),STE),(RTB(4),CTS),(RTB(5),STS)

      LOGICAL SHARPE, SHARPS
 
      XLIM =A(NOEL,10)
      NP = A(NOEL,11)
      BBM(1)  =A(NOEL,12)*SCAL

      XE = A(NOEL,20)
      DLE(1) = A(NOEL,21)

      XLS = A(NOEL,40)
      DLS(1) = A(NOEL,41)

      DO 8 I=1,6
        CE(I) =A(NOEL,30+I)
        CS(I) =A(NOEL,50+I)
 8    CONTINUE
 
C----------- Champ DE FUITE
      SHARPE=DLE(1) .EQ. 0.D0
      SHARPS=DLS(1) .EQ. 0.D0
      IF(SHARPE) XE=0.D0
      IF(SHARPS) XLS=0.D0
      XI = ZERO
      XL=XLIM
      XLIM = XL + XE + XLS
      XF = XLIM
      XS = XL + XE
      CTE=UN
      STE=ZERO
      CTS=UN
      STS=ZERO
 
C----- SHARP EDGE => INTEGR STOPPE SUR DR. DE COUPURE
      IF(SHARPE) THEN
C------- Correction for entrance wedge
          CALL INTEG1(ZERO,ZERO)
C------------------------------------
        IDRT = -1
        CA(1)=CTE
        SA(1)=STE
        CM(1)=-XE*CA(1)
      ELSE
        DE(1,1)=-BBM(1)/DLE(1)
        DO 10 I=2,4
          DE(1,I)=-DE(1,I-1)/DLE(1)
 10     CONTINUE
      ENDIF

      IF(SHARPS) THEN
C------- Correction for exit wedge
          CALL INTEG2(ZERO,ZERO)
C------------------------------------
        IF(IDRT .EQ. -1) THEN
          IDRT = 2
        ELSE
          IDRT = 1
        ENDIF
        CA(2)=UN
        SA(2)=ZERO
        CM(2)=-XS*CA(2)
      ELSE
        DS(1,1)=-BBM(1)/DLS(1)
        DO 11 I=2,4
          DS(1,I)=-DS(1,I-1)/DLS(1)
 11     CONTINUE
      ENDIF
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,100) ' UNDULATOR',XL,NP
 100    FORMAT(/,5X,' +++++  ',A10,'  : '
     >        ,//,15X,' Length              = ',F10.3,'  cm'
     >         ,/,15X,' Number  of  periods = ',I8)
        WRITE(NRES,103) ' MAX     ', BBM(1)
 103    FORMAT(15X,' B-',A,'  =',F12.6,'  KG ')
 
        WRITE(NRES,104) 'D''ENTREE'
 104    FORMAT(/,15X,' FACE  ',A)
        WRITE(NRES,101) XE,DLE(1)
 101    FORMAT(15X,' DX = ',F10.3,'    LAMBDA = ',F10.3)
        IF( .NOT. SHARPE ) WRITE(NRES,132) (CE(I),I=1,6)
 132    FORMAT(15X,' Fringe  field  coefficients :',/,16X,6F9.5)
C 132    FORMAT(15X,' COEFFICIENTS DE Champ DE FUITE:',/,16X,6F9.5)
 
        WRITE(NRES,104) 'DE  SORTIE'
        WRITE(NRES,101) XLS,DLS(1)
        IF( .NOT. SHARPS ) WRITE(NRES,132) (CS(I),I=1,6)

        IF( SHARPE .OR. SHARPS ) 
     >     WRITE(NRES,FMT='(/,''  ***  Warning : sharp edge '',
     >     ''model entails vertical wedge focusing simulated with'',
     >     /,17X,'' first order kick  ***'')')

      ENDIF

      CALL UNDUL2(XL,NP)

      RETURN
      END
