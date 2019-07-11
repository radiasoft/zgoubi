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
      SUBROUTINE EL2TUB(SCAL,XL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.AIM.H"     ! COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI

      EQUIVALENCE (EN,D), (EB1,V21), (EB2,OM), (EG1,X0)

      X1 =A(NOEL,10)
      D =A(NOEL,11)
      X2 =A(NOEL,12)
      RO =A(NOEL,13)
      V1 =A(NOEL,20)*SCAL
      V2 =A(NOEL,21)*SCAL
      XLIM=X1+D+X2

      IF(NRES.GT.0) THEN
        WRITE(NRES,110) XLIM,D,RO,V1,V2
 110    FORMAT(/,25X,' --- LENTILLE ELECTROSTATIQUE A 2 TUBES ---'
     >        ,/,30X,' Length                    = ',1P,G12.5,' m'
     >        ,/,30X,' DISTANCE  between  tubes  = ',G12.5,' m'
     >        ,/,30X,' Radius                    = ',G12.5,' m'
     >        ,/,30X,' Potential  1              = ',G12.5,' V'
     >        ,/,30X,' Potential  2              = ',G12.5,' V')

C------- DISTANCES FOCALES
        IF(V1*V2 .NE. 0.D0) THEN
          OM=1.318D0/RO
          GA=V1/V2
          G2=SQRT(GA)
          G4=SQRT(G2)
          ZM=LOG(GA)/4.D0/OM
          EZ=EXP(OM*ZM)
          TM=-2.D0*OM*(EZ-1.D0/EZ)/(EZ+1.D0/EZ)
          AA=-((1.D0+G2)/(1.D0-G2))**2*(2.D0+(1.D0+GA)/
     >          (1.D0-GA)*LOG(GA))/PI/OM
          AK=SQRT(1.D0+3.D0/16.D0*TM*TM*AA*AA)
          SF=SIN(PI*AK)
          F2=AA*AK/SF/G4
          F1=F2*G4*G4
          WRITE(NRES,111) F1,F2
 111      FORMAT(/,29X,' DISTANCES  FOCALES:'
     >        ,/,30X,' F1 = ',G12.4,' M,  F2 = ',G12.4,' M')
        ENDIF
      ENDIF

C----- change unites: XLIM:cm, R:cm, V1-2:MeV, D:cm
      XLIM=XLIM*100.D0
      RO=RO*100.D0
      X0=(X1+0.5D0*D)*100.D0

C----- Omega=1.318/RO
      OM=1.318D0/RO
      V1=V1*1.D-6
      V2=V2*1.D-6
      D=D*100.D0

      IF(V1 .EQ. V2) KFLD=0
      V21=V2-V1
      XI = 0.D0
      XF=XLIM
      XL=XLIM

      RETURN
      END
