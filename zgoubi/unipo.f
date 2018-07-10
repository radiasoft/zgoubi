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
      SUBROUTINE UNIPO(XX,Y,Z,BRI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.AIM.H"     ! COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CHAVE_2.H"     ! COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      INCLUDE "C.DDEXYZ.H"     ! COMMON/DDEXYZ/ DE(3,3),DDE(3,3,3)
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
 
C----- D=DISTANCE TUBES, VDR=(V2-V1)/(D.RO), OR=Omega/RO, X0=refX
      EQUIVALENCE (EN,D), (EB1,VDR),(EB2,OR),(EG1,X0),(EG2,X22)
 
      PARAMETER (MDX=6)
      DIMENSION EX(MDX)
      DIMENSION ER(2),DER(2,2),DDER(2,2,2)
 
      CALL EAX3TU(OR,XX-X0,VDR*BRI,D,X22,EX) 
      R2  =Y*Y + Z*Z
      R   =SQRT(R2)
      CALL BAXBXR(EX,R,R2,
     >                    ER,DER,DDER)
      CALL BXRXYZ(ER,DER,DDER,Y,Z,R,2,
     >                                   E,DE,DDE)
 
      RETURN
      END
