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
C  USA
C  -------
      SUBROUTINE SOLENO(KUASEX,LMNT,SCAL,
     >                                   XL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) LMNT(*)
      INCLUDE "C.AIM.H"     ! COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      INCLUDE "C.CHAFUI.H"     ! COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.DROITE.H"     ! COMMON/DROITE/ CA(9),SA(9),CM(9),IDRT
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      PARAMETER(MPOL=10)
      INCLUDE "C.MULTPE.H"     ! COMMON/MULTPE/ EM(MPOL),QLE(MPOL),QLS(MPOL)
      INCLUDE "C.MULTPL.H"     ! COMMON/MULTPL/ BM(MPOL),DLE(MPOL),DLS(MPOL),DE(MPOL,MCOEF),DS(MPOL,MCOEF),RTB(MPOL)
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE "MAXTRA.H"
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM

      DIMENSION  AREG(2),BREG(2),CREG(2)
      CHARACTER(20) TYP(2)
      DATA TYP / 'axial field model', 'elliptic integrals' /

      XL =A(NOEL,10)
      RO =A(NOEL,11)
      BO =A(NOEL,12)*SCAL
C MODL=1 (default) for axial model, MODL=2 for E,K,PI integral model
      MODL = NINT(A(NOEL,13))
      XE = A(NOEL,20)
      XLS= A(NOEL,21)

      RO2 = RO*RO
C      BOSQ = BO * SQRT(XL*XL+4.D0*RO2)/(2.D0*XL)
C      CALL SOLEN2(MODL,BOSQ,RO2)
C FM - July 2015
      BO2 = BO /2.D0
      CALL SOLEN2(MODL,BO2,RO2)

      IF(NRES.GT.0) THEN
        WRITE(NRES,100) LMNT(KUASEX),XL,RO
 100    FORMAT(/,5X,' -----  ',A10,'  : ', 1P
     >  ,/,15X,' Length  of  element  : ',G12.4,'  cm'
     >  ,/,15X,' Inner radius  RO =',G12.4,'  cm')
        WRITE(NRES,103) 'CNTRL',BO,' kG ;  ',BO/BORO
     >  ,BO*XL/(2.D0*BORO)
 103    FORMAT(15X,' B-',A,'  =',1P,G12.4,1X,A,' K=B/BORO = '
     >  ,G12.4,' /m ;   theor. angle  BL/(2*BORO) = ',G12.4,' rad')
        WRITE(NRES,145) XE,XLS
 145    FORMAT(15X,' Entrance and exit integration extents : '
     >  /,20X,'   XE =',G12.4,' cm,   XS =',1P,G12.4,' cm')
        WRITE(NRES,FMT='(15X,'' MODL = '',I2,
     >  '' -> Solenoid model is '',A)') MODL,TYP(MODL)
        IF(KSPN.EQ.1) THEN
          WRITE(NRES,FMT='(15X,'' Spin X-rotation angle phi=(1+G)*B0.L'
     >   //'/Brho ) = '',F10.4,''.  Strength phi/180 = '',F10.4,''%'')')
     >   (1.D0+G)*BO*XL/BORO*DEG,(1.D0+G)*BO*XL/BORO/PI*1.D2
        ENDIF
      ENDIF

      XI = 0.D0
      XLIM = XL + XE + XLS
      XF = XLIM
      XS = XLIM - XLS
      IF(BO .EQ. 0.D0) KFLD=0

      CALL CHXC1R(
     >            KPAS)
      IF(KPAS.GE.1) THEN
        IF(XE+XLS .GE.XL) THEN
          CALL ENDJOB(' Pgm soleno. Entrance/body/exit step size mode '
     >    //'is not compatible with fringe fields overlapping in body.'
     >    //' Instead, use explicit single step size value. ',-99)
        ENDIF
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
