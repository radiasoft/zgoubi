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
      SUBROUTINE OPTIMP(LUN,NOEL,F0,PHY,PHZ,AKL,CSTRN,RPRM,R,
     >                                                     PP0) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,6), F0(6,6),AKL(3)

      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE 'MXLD.H'
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LABEL
      INCLUDE "C.LABEL.H"     ! COMMON/LABEL/ LABEL(MXL,2)
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI

      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLE

      CHARACTER(LBLSIZ) LBL1, LBL2
      LOGICAL EMPTY

      PARAMETER (PI = 4.D0*ATAN(1.D0))
C      PARAMETER (PI = 3.1415926535898d0)

      CALL SCUMR(
     >            XL,SCUM,TCUM)

      CALL ZGKLEY( 
     >            KLE)
      LBL1 = LABEL(NOEL,1)
      IF(EMPTY(LBL1)) LBL1 = '-'
      LBL2 = LABEL(NOEL,2)
      IF(EMPTY(LBL2)) LBL2 = '-'

      CALL ZGIPAS( 
     >            IPASS,NRBLT)

C      WRITE(LUN,104) -F0(1,2), F0(1,1), -F0(3,4), F0(3,3) 
C FM Mar 2015
C                     +alfY     betY    +alfZ    betZ
      WRITE(LUN,104) F0(1,2), F0(1,1), F0(3,4), F0(3,3) 
     >, F0(5,6), F0(5,5)
     >, F0(1,6), F0(2,6), F0(3,6), F0(4,6)
     >, PHY/(2.D0*PI), PHZ/(2.D0*PI), SCUM, NOEL
C     >, PHY, PHZ, SCUM, NOEL
     >, F(2,1)*0.01D0, F(3,1), F(4,1)*0.01D0, F(5,1)
     >, kle, LBL1, LBL2, F(6,1)*0.01D0, (AKL(I), I=1,3)
     >, CSTRN, RPRM, '   ! optimp.f', IPASS,DPREF
     >, R(1,1),R(1,2),R(2,1),R(2,2)
     >, R(3,3),R(3,4),R(4,3),R(4,4)
     >, R(5,1),R(5,2),R(5,3),R(5,4),R(5,6)
 104  FORMAT(1P,13(E13.5,1X),1X,I5,4(1X,E13.5),3(1X,A)
     >,6(1X,E13.5),A,1X,I0,1X,14(E13.5,1X))
C FM 14 Apr. 14
C     >,5(1X,E13.5),A,/)

      PP0 = F(1,1)

      RETURN
      END
