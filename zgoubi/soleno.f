C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <meot@lpsc.in2p3.fr>
C  Service Acc�lerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE SOLENO(KUASEX,LMNT,SCAL,
     >                                   XL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER LMNT(*)*(*)
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
     >,QE(MPOL,6),QS(MPOL,6),RTQ(MPOL)
      COMMON/MULTPL/ BM(MPOL),DLE(MPOL),DLS(MPOL)
     >,DE(MPOL,10),DS(MPOL,10),RTB(MPOL)
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
 
      DIMENSION  AREG(2),BREG(2),CREG(2)
 
          XL =A(NOEL,10)
          RO =A(NOEL,11)
          BO =A(NOEL,12)*SCAL
          XE = A(NOEL,20)
          XLS= A(NOEL,21)
 
          IF(NRES.GT.0) THEN
            WRITE(NRES,100) LMNT(KUASEX),XL,RO
 100        FORMAT(/,5X,' -----  ',A10,'  : ', 1P
     >      ,/,15X,' Length  of  element  : ',G12.4,'  cm'
     >      ,/,15X,' Inner radius  RO =',G12.4,'  cm')
            WRITE(NRES,103) 'CNTRL',BO,' kG'
 103        FORMAT(15X,' B-',A,'  =',1P,G12.4,1X,A6)
            WRITE(NRES,145) XE,XLS
 145        FORMAT(20X,' XE =',G12.4,' cm,   XS =',1P,G12.4,' cm')
          ENDIF

          XI = 0.D0
          XLIM = XL + XE + XLS
          XF = XLIM
          XS = XLIM - XLS
          IF(BO .EQ. 0.D0) KFLD=0

      CALL CHXC1R(
     >            KPAS)
      IF(KPAS.GE.1) THEN
        AREG(1)=1.D0
        BREG(1)=0.D0
        CREG(1)=-2.D0*XE
        AREG(2)=1.D0
        BREG(2)=0.D0
        CREG(2)=-2.D0*XS+XLIM
        CALL INTEG6(AREG,BREG,CREG)
      ENDIF

      RETURN
      END
