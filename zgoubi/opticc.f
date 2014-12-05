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
      SUBROUTINE OPTICC(NOEL,KOPIMP,PRDIC,OKCPLD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,6), F0(6,6), AKL(3)
      logical OKCPLD

      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLE
      LOGICAL PRDIC
      LOGICAL OKLNO, OKLNOI
      SAVE OKLNO
      SAVE LNOPTI

      IORD = 1
      IFOC = 0
      KWR = 0
      CALL MATRIC(IORD,IFOC,KWR,OKCPLD)
      CALL MATRI1(
     >            R)

      CALL BEAMAT(R,PRDIC,OKCPLD,
     >                           F0,PHY,PHZ,Cstrn)
      CALL BEAIMP(F0,PHY,PHZ)

      CALL ZGKLEY(KLE)
      IF(KLE.EQ.'AGSMM') THEN
         CALL AGSMKL(
     >        AL, AK1, AK2, AK3)
         AKL(1) = AK1 * AL *1.D-2
         AKL(2) = AK2 * AL *1.D-2
         AKL(3) = AK3 * AL *1.D-2
      ELSEIF(KLE.EQ.'AGSQUAD') THEN
         CALL AGSQKL(
     >        AL, AK1)
         AKL(1) = 0.D0
         AKL(2) = AK1 * AL *1.D-1
      ELSEIF(KLE.EQ.'MULTIPOL') THEN
         CALL MULTKL(
     >        AL, AK1, AK2, AK3)
         AKL(1) = AK1 * AL *1.D-2
         AKL(2) = AK2 * AL *1.D-2
         AKL(3) = AK3 * AL *1.D-2
      ELSE
         AKL = 0.D0
      ENDIF

c Print to zgoubi.OPTICS.out
c      IF(KOPIMP.EQ.1) THEN

      IF(OKLNO) 
     > CALL OPTIMP(LNOPTI,NOEL,F0,PHY,PHZ,AKL,Cstrn)      

c Store for later print by TWISS to zgoubi.TWISS.out
c      ELSEIF(KOPIMP.EQ.2) THEN

c      ENDIF

      RETURN

      ENTRY OPTIC2(
     >             OKLNOI,LNOPTJ)
      OKLNO = OKLNOI
      LNOPTI = LNOPTJ
      RETURN
      END
