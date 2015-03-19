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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE COILS(SCAL,
     >                              XLT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.AIM.H"     ! COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      INCLUDE "C.CHAFUI.H"     ! COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
C      COMMON/CHAFUI/ XE,XS,CE(6),CS(6),QCE(6),QCS(6)
      PARAMETER (MXCOIL=8)
      INCLUDE "C.COIL.H"     ! COMMON/COIL/ XLI(MXCOIL),RI(MXCOIL),BI(MXCOIL),DIST(MXCOIL),MS
      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.DROITE.H"     ! COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
 
      DIMENSION  AREG(2),BREG(2),CREG(2)
 
      MS = NINT(A(NOEL,2))
      XLT = 0.D0
      BO = 0.D0
       
      M=0
      M1 = M+1
      DIST(M1) = 0.D0
      XLI(M1) = A(NOEL,3+4*M)
      RI(M1) =  A(NOEL,4+4*M)
      BI(M1) =  A(NOEL,5+4*M)*SCAL
      DO 1 M = 1, MS-1         
          M1 = M+1
          IF(M1.GT.MXCOIL) 
     >      CALL ENDJOB(' Too  many  coils, max is ',MXCOIL)
          DIST(M1) =A(NOEL,6+4*(M-1))
          XLI(M1) = A(NOEL,3+4*M)
          RI(M1) =  A(NOEL,4+4*M)
          BI(M1) =  A(NOEL,5+4*M)*SCAL
          XLT = XLT + DIST(M)
          BO = BO + BI(M) 
 1    CONTINUE
      XLT = XLT + DIST(MS)
      XLT = XLT + XLI(1)/2.D0 + XLI(MS)/2.D0 

      XE = A(NOEL,7+4*MS)
      XLS = A(NOEL,8+4*MS)

        IF(NRES.GT.0) THEN
          WRITE(NRES,*) 'NUMBER  OF  COILS :',MS
          DO 2 M = 1,MS 
            WRITE(NRES,100) 'Coil  #',M,DIST(M),XLI(M),RI(M),BI(M)
 100        FORMAT(/,5X,A,I2,1P, 
     >      /,15X,' Distance  from  previous  coil  = ',G12.4,'  cm'
     >      /,15X,' Length  of  element  = ',G12.4,'  cm'
     >      /,15X,' Radius   RO          = ',G12.4,'  cm'
     >      /,15X,' Asymptotic field     = ',G12.4,'  kG')
 2        CONTINUE

          WRITE(NRES,145) XE,XLS
 145      FORMAT(/,20X,' XE =',F7.3,'  cm,  XS =',F7.3,'  cm')
        ENDIF

        XI = 0.D0
        XLIM = XLT + XE + XLS
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
