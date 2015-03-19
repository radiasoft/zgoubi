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
      SUBROUTINE ELCMII(SCAL,
     >                       XL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.AIM_2.H"     ! COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      INCLUDE "C.CHAFUI.H"     ! COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      PARAMETER (MDR=9)
      INCLUDE "C.DROITE.H"     ! COMMON/DROITE/ CA(MDR),SA(MDR),CM(MDR),IDRT
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
  
      R1 = A(NOEL,10)
      R2 = A(NOEL,11)
      AT = A(NOEL,12)
      D = A(NOEL,13)

      V21 =A(NOEL,20)*SCAL
      V32 =A(NOEL,21)*SCAL

C------ NORMALLY RB-R2 AND R1-RA ~ 3D
      
      IF(NRES.GT.0) THEN
        WRITE(NRES,110) R1,R2,AT,D,V21,V32
 110    FORMAT(/,40X,' --- TRANSAXIAL  MIRROR ---',/,1P
     >        ,/,30X,' Radius AT SLIT 1              (R1) = ',G12.5,' M'
     >        ,/,30X,' Radius AT SLIT 2              (R2) = ',G12.5,' M'
     >        ,/,30X,' TOTAL ANGLE         (AT) = ',G12.5,' RAD'
     >        ,/,30X,' GAP                 (D)  = ',G12.5,' M'
     >        ,/,30X,' V2 - V1                  = ',G12.5,' V'
     >        ,/,30X,' V3 - V2                  = ',G12.5,' V')
      ENDIF
 
C----- CHANGE UNITS: R1,R2,D -> CM 
      R1 = R1 * 100.D0
      R2 = R2 * 100.D0
      D = D * 100.D0
C----- V1-3 -> V/D, MEV/CM
      V21D = V21*1.D-6/D
      V32D = V32*1.D-6/D

C------ RM IS IN CM
      RM = 0.5D0*(R1+R2) 
      AE = 0.D0
      AS = 0.D0
C------ XI, XF IS IN RAD
      XI = -AT/2.D0
      XF = AT/2.D0

      IF(V21*V21 + V32*V32 .EQ. 0.D0) KFLD=0

      CALL ELCMIW(D,R1,R2,V21D,V32D)

      XL = 0.0

      RETURN
      END
