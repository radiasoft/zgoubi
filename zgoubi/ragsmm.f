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
C  François Meot <fmeot@bnl.gov>
C  Brookhaven National Laboratory
C  C-AD, Bldg 911
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE RAGSMM(NDAT,NOEL,MXL,A,ND)
C     -------------------------
C     READS DATA FOR AGS DIPOLE
C     -------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MXL,*)
      CHARACTER(130) TXT
      PARAMETER (MSS=6)
      CHARACTER(80) STRA(MSS)
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLE

C --- IL
      LINE = 1
      READ(NDAT,*,ERR=90,END=90) A(NOEL,1)

C mod = option flag (for K1, K2 method)
C                  MOD, dL/L, gap, db0/b0, db1/b1, db2/b2
      LINE = LINE + 1
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,10+I-1),I=1,6)
C # of backleg windings, # of turns and I for each winding
      LINE = LINE + 1
      READ(NDAT,*,ERR=90,END=90) ABLW,(A(NOEL,20+I),I=1,2*INT(ABLW))
      IF(INT(ABLW) .GT. 2) STOP ' SBR RAGSMM. # of blwg cannot exceed 2'
      A(NOEL,20) = ABLW
C----- CHP FUITE ENTREE
C        XE, LambdaE, ff2,ff3
      LINE = LINE + 1
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=30,33)
      LINE = LINE + 1
      READ(NDAT,*,ERR=90,END=90) II,(A(NOEL,40+I),I=1,II)
      A(NOEL,40) = II
C----- CHP FUITE SORTIE
      LINE = LINE + 1
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,I),I=50,53)
      LINE = LINE + 1
      READ(NDAT,*,ERR=90,END=90) II,(A(NOEL,60+I),I=1,II)
      A(NOEL,60) = II
C----- Rotation of multipole components
      LINE = LINE + 1
      READ(NDAT,*,ERR=90,END=90) (A(NOEL,70+I-1),I=1,3)

C----- Step
      ND = 80
      CALL STPSIZ(NDAT,NOEL,ND,
     >                         A)

C KPOS=3 : XCE, YCE, ALE
C KPOS=4 : X-shft, Y-shft, Z-rot, Z-shft, Y-rot
      LINE = LINE + 1
      READ(NDAT,FMT='(A)',ERR=90,END=90) TXT
      CALL STRGET(TXT,MSS,
     >                    NS,STRA)
      READ(STRA(1),*,ERR=66) II
      A(NOEL,90) = II
      DO KK = 2, NS
        READ(STRA(KK),*,ERR=66) TEMP
        A(NOEL,91+KK-2) = TEMP
      ENDDO

      RETURN

 66   CONTINUE
      RETURN

 90   CONTINUE
      CALL ZGKLEY(
     >            KLE)
      CALL ENDJOB('*** Pgm ragsmm, keyword '//KLE//' : '//
     >'input data error, at line #',LINE)
      RETURN
      END
