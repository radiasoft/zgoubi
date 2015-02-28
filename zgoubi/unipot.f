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
      SUBROUTINE UNIPOT(SCAL,XL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.AIM.H"     ! COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
 
      EQUIVALENCE (EN,D), (EB1,VDR),(EB2,OR),(EG1,X0),(EG2,X22)
 
      X1 =A(NOEL,10)
      D =A(NOEL,11)
      X2 =A(NOEL,12)
      X3 =A(NOEL,13)
      RO =A(NOEL,14)
      V1 =A(NOEL,20) *SCAL
      V2 =A(NOEL,21) *SCAL
      XL=X1+D+X2+D+X3

      IF(NRES.GT.0) THEN
        WRITE(NRES,110) XL,D,RO,V1,V2
 110    FORMAT(/,25X,' --- UNIPOTENTIAL LENSE ---',1P
     >        ,/,30X,' Length                    = ',G12.5,' m'
     >        ,/,30X,' Distance  between  tubes  = ',G12.5,' m'
     >        ,/,30X,' Radius                    = ',G12.5,' m'
     >        ,/,30X,' Potential  1              = ',G12.5,' V'
     >        ,/,30X,' Potential  2              = ',G12.5,' V')
      ENDIF
 
C----- change unites: XL, R:cm, V1-2:MeV, D:cm
      XL=XL*100.D0
      RO=RO*100.D0
      X22=.5D0*X2*100.D0
      D=D*100.D0
      X0=X1*100.D0+ D+X22
      V1=V1*1.D-6
      V2=V2*1.D-6

C----- Omega=1.318, OR=Omega/RO
      OR=1.318D0/RO
 
      IF(V1 .EQ. V2) KFLD=0
      VDR=(V2-V1)/D/RO
      XI = 0.D0
      XLIM=XL
      XF=XLIM
                           
      RETURN
      END
