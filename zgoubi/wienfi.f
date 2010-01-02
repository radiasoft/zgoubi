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
      SUBROUTINE WIENFI(SCAL,XL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
C      COMMON/CHAFUI/ XE,XS,CE(6),CS(6),QCE(6),QCS(6)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QEL,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      COMMON/DROITE/ IDRT,CA(9),SA(9),CM(9)
      COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      PARAMETER(MPOL=10)
      COMMON/MULTPE/ EM(MPOL),QLE(MPOL),QLS(MPOL)
     >,QE(MPOL,MCOEF),QS(MPOL,MCOEF),RTQ(MPOL)
      COMMON/MULTPL/ BM(MPOL),DLE(MPOL),DLS(MPOL)
     >,DE(MPOL,MCOEF),DS(MPOL,MCOEF),RTB(MPOL)
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
 
      EQUIVALENCE (RTB(1),CTE),(RTB(2),STE),(RTB(4),CTS),(RTB(5),STS)
      EQUIVALENCE
     >  (RTQ(1),CTEQ),(RTQ(2),STEQ),(RTQ(4),CTSQ),(RTQ(5),STSQ)
      LOGICAL SHARPE, SHARPS
 
      CHARACTER*10 TYP(2)
      DATA TYP / 'HORIZONTAL',' VERTICAL ' /
 
      XL =A(NOEL,10)
      E  =A(NOEL,11)*SCAL
      B  =A(NOEL,12)*SCAL
      KHV =A(NOEL,13)
      XE = A(NOEL,20)
      QLE(1) = A(NOEL,21)
      DLE(1) = A(NOEL,22)
      XLS = A(NOEL,50)
      QLS(1) = A(NOEL,51)
      DLS(1) = A(NOEL,52)
      DO 8 I=0,5
        I1=I+1
        QCE(I1) =A(NOEL,30+I)
        CE(I1) =A(NOEL,40+I)
        QCS(I1) =A(NOEL,60+I)
        CS(I1) =A(NOEL,70+I)
 8    CONTINUE
 
      TE=0.D0
      TS=0.D0
      CTE=COS(TE)
      STE=SIN(TE)
      CTS=COS(TS)
      STS=SIN(TS)
      TEQ=0.D0
      TSQ=0.D0
      CTEQ=COS(TEQ)
      STEQ=SIN(TEQ)
      CTSQ=COS(TSQ)
      STSQ=SIN(TSQ)
 
C----------- Champ DE FUITE
      SHARPE=QLE(1)*DLE(1) .LE. 0.D0
      SHARPS=QLS(1)*DLS(1) .LE. 0.D0
      IF(SHARPE) THEN
C------- Correction for entrance wedge
        FINTE = XE
        XE=0D0
        GAPE = -DLE(1)
C        GAPE = -AMAX1(QLE(1),DLE(1))
        IF(NRES.GT.0) 
     >  WRITE(NRES,FMT='(''Entrance hard edge is to be implemented'')')
        CALL INTEG1(-TE,FINTE,GAPE)
C        QLE(1)=0.D0
C        DLE(1)=0.D0
      ENDIF
      IF(SHARPS) THEN
C------- Correction for exit wedge
        FINTS = XLS
        XLS=0D0
C        GAPS = -AMAX1(QLS(1),DLS(1))
        GAPS = -DLS(1)
        IF(NRES.GT.0) 
     >  WRITE(NRES,FMT='(''Exit hard edge is to be implemented'')')
        CALL INTEG2(-TS,FINTS,GAPS)
C        QLS(1)=0.D0
C        DLS(1)=0.D0
      ENDIF
 
C---------------------------------------------------
C        KHV = 1 (H SEPARATION) OR 2 ( V SEPARATION)
C        E  ( V/M )  = E-FIELD
C        B  ( T )    = B-FIELD
C---------------------------------------------------
 
C      P0 = BORO*CL*1.D-9*Q/QEL
C      P0 = BORO*CL*1.D-9*Q
      P0 = BORO*CL9*Q
      IF(NRES.GT.0) THEN
        WRITE(NRES,110) TYP(KHV),XL,E,B,AM,Q,P0/SQRT(P0*P0+AM*AM)
 110    FORMAT(/,25X,' ------ SEPARATEUR ELECTROSTATIQUE ------'
     >        ,/,40X,A
     >        ,/,30X,' Length                  = ',1P,G12.5,' m'
     >        ,/,30X,' E                       = ',G12.5,' V/m'
     >        ,/,30X,' B                       = ',G12.5,' T'
     >      ,/,/,30X,' Mass  of  particles     = ',G12.5,' MeV/c2'
     >      ,/,/,30X,' Charge of  particles    = ',G12.5,' C'
     >        ,/,30X,' Reference  beta         = ',G12.5)
        IF(B .NE. 0.D0) WRITE(NRES,111) E/(-B*CL)
 111    FORMAT(30X,' BETA  WANTED = -E/c.B =',1P,G12.5)
 
        WRITE(NRES,104) 'D''ENTREE'
 104    FORMAT(/,15X,' FACE  ',A)
        WRITE(NRES,101) XE,QLE(1),DLE(1)
 101    FORMAT(15X,' DX = ',F7.3,',  LAMBDA-E, B = ',2F7.3,' CM')
        IF(QLE(1) .NE. 0.D0) WRITE(NRES,132) 'E',(QCE(I),I=1,6)
 132    FORMAT(15X,' COEFFICIENTS DE ',A,' :',6F9.5)
        IF(DLE(1) .NE. 0.D0) WRITE(NRES,132) 'B',(CE(I),I=1,6)
 
        WRITE(NRES,104) 'DE  SORTIE'
        WRITE(NRES,101) XLS,QLS(1),DLS(1)
        IF(QLS(1) .NE. 0.D0) WRITE(NRES,132) 'E',(QCS(I),I=1,6)
        IF(DLS(1) .NE. 0.D0) WRITE(NRES,132) 'B',(CS(I),I=1,6)
 
      ENDIF
 
C----- change unites: XL->cm, B->kG, E->MeV/cm
      XL=XL*100.D0
      XI = 0.D0
      XLIM = XL + XE + XLS
      XF = XLIM
      XS = XL + XE
 
      EM(1)=E*1.D-8
      BM(1)=B*10.D0
      IF    (KHV .EQ. 1) THEN
       BM(6)=0.D0
       EM(6)=.5D0*PI
      ELSEIF(KHV .EQ. 2) THEN
       EM(6)=0.D0
       BM(6)=-.5D0*PI
      ENDIF
 
      IF(B .EQ. 0.D0) KFLD=LC
      IF(E .EQ. 0.D0) KFLD=MG
      IF(E*E+B*B .EQ. 0.D0) KFLD=0
 
C----- DEFINITION DES ANGLES DE COIN
C      (=0 DEG. DANS CETTE VERSION)
      IF(SHARPE) THEN
        IDRT = -1
        CA(1)=1.D0
        SA(1)=0.D0
        CM(1)=-XE
      ELSE
        QE(1,1)=-EM(1)/QLE(1)
        DE(1,1)=-BM(1)/DLE(1)
        DO 10 I=2,3
          QE(1,I)=-QE(1,I-1)/QLE(1)
          DE(1,I)=-DE(1,I-1)/DLE(1)
 10     CONTINUE
      ENDIF
      IF(SHARPS) THEN
        IF(IDRT .EQ. -1) THEN
          IDRT = 2
        ELSE
          IDRT = 1
        ENDIF
        CA(2)=1.D0
        SA(2)=0.D0
        CM(2)=-XS
      ELSE
        QS(1,1)=-EM(1)/QLS(1)
        DS(1,1)=-BM(1)/DLS(1)
        DO 11 I=2,3
          QS(1,I)=-QS(1,I-1)/QLS(1)
          DS(1,I)=-DS(1,I-1)/DLS(1)
 11     CONTINUE
      ENDIF

      IF( SHARPE .OR. SHARPS ) 
     >     WRITE(NRES,FMT='(/,''  ***  Warning : sharp edge '',
     >    ''model entails vertical wedge focusing simulated with'',
     >    /,17X,'' first order kick  ***'')')

      RETURN
      END
