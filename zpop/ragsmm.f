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
      SUBROUTINE RAGSMM(NDAT,NOEL,MXL,A,ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MXL,*)
      CHARACTER(130) TXT
C     -------------------------
C     READS DATA FOR AGS DIPOLE
C     -------------------------
C                  IL
      LINE = 1
      READ(NDAT,*) A(NOEL,1)

c mod = option flag (for K1, K2 method)
C                  MOD, dL/L, gap, db0/b0, db1/b1, db2/b2
      LINE = 2
      READ(NDAT,*) (A(NOEL,10+I-1),I=1,6)
C       # of turns and I for each winding   
      LINE = 3
      READ(NDAT,*) NBLWG,(A(NOEL,20+I),I=1,2*NBLWG)
      IF(NBLWG.GT.2) STOP ' SBR ragsmm. # of blwg cannot exceed 2'
      A(NOEL,20) = NBLWG
C----- CHP FUITE ENTREE
C        XE, LambdaE, ff2,ff3
      LINE = 4
      READ(NDAT,*) (A(NOEL,I),I=30,33)
      LINE = 5
      READ(NDAT,*) II,(A(NOEL,40+I),I=1,II)
      A(NOEL,40) = II
C----- CHP FUITE SORTIE
      LINE = 6
      READ(NDAT,*) (A(NOEL,I),I=50,53)
      LINE = 7
      READ(NDAT,*) II,(A(NOEL,60+I),I=1,II)
      A(NOEL,60) = II
C----- Rotation of multipole components
      LINE = 8
      READ(NDAT,*) (A(NOEL,70+I-1),I=1,3)

      LINE = 10
C----- KPOS, XCE, YCE, ALE
      READ(NDAT,*,err=99) II,(A(NOEL,I),I=91,93)

 99   RETURN
      END
