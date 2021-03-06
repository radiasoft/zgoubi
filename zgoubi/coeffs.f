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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE COEFFS(IOPT,IORD,R,T,IREF,
     >                                     F0,Cstrn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*),T(6,6,*),F0(6,*)
C     ----------------------------------------------------------------------
C     Transfer  coefficients
C     Called by ff during FIT process. IOPT=1 if fitting beam matrix coeffs.
C     ----------------------------------------------------------------------
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)

      DIMENSION TX3(5,6) , TX4(5,6)
      LOGICAL PRDIC
      DATA PRDIC / .FALSE. /
C      LOGICAL FITING

      IF    (IORD .EQ. 1) THEN
        IT1 = 1 + 11 * (IREF-1)
        IT2 = IT1+3
        IT3 = IT1+4
        CALL REFER(1,1,0,IT1,IT2,IT3)
        CALL MAT1(IT1,F,IMAX,
     >                R,T)
        CALL REFER(2,1,0,IT1,IT2,IT3)
      ELSEIF(IORD .EQ. 2) THEN
        CALL REFER(1,2,0,1,6,7)
        CALL MAT2(
     >            R,T,TX3,TX4)
        CALL REFER(2,2,0,1,6,7)
      ENDIF

      CALL MKSA(IORD,R,T,TX3,TX4)

      PRDIC = .FALSE.
      IF(IOPT.EQ.1) CALL BEAMAT(R,PRDIC,.FALSE.,
     >                                          F0,PHY,PHZ,CSTRN,RPRM)

      RETURN
      END
