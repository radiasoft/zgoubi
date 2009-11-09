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
      SUBROUTINE UNDULI(SCAL,
     >                       XL)
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
     >,DE(MPOL,MCOEF),DS(MPOL,MCOEF),RTB(MPOL)
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
 
      EQUIVALENCE (RTB(1),CTE),(RTB(2),STE),(RTB(4),CTS),(RTB(5),STS)

C      LOGICAL SHARPE, SHARPS
 
      XL =A(NOEL,10)
      BBM(1)  =A(NOEL,11)*SCAL
      AA = A(NOEL,12)
      BB = A(NOEL,13)

      XLIM=XL
      XI = -XLIM/2.D0
      XF = XLIM/2.D0
      CTE=UN
      STE=ZERO
      CTS=UN
      STS=ZERO
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,100) ' UNDULATOR :  B(x)=Bo.sin(x/B)/(1+(x/A)^8))',XL
 100    FORMAT(/,5X,' -----  ',A10,'  : ', 1P, 
     >        //,15X,' Length              = ',G12.4,'  cm')
        WRITE(NRES,103) ' MAX     ', BBM(1), AA, BB
 103    FORMAT(15X,1P,' B-',A,'  =',G12.4,'  kG ', 
     >        /15X   ,' A        =',G12.4,'  cm',
     >        /15X   ,' B        =',G12.4,'  cm')
      ENDIF
 
      CALL UNDUL1(AA,BB)

      RETURN
      END
