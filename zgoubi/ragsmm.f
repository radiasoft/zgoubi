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
      PARAMETER (MSS=6)
      CHARACTER(80) STRA(MSS)
C     -------------------------
C     READS DATA FOR AGS DIPOLE
C     -------------------------
C                  IL
      READ(NDAT,*) A(NOEL,1)

C mod = option flag (for K1, K2 method)
C                  MOD, dL/L, gap, db0/b0, db1/b1, db2/b2
      READ(NDAT,*) (A(NOEL,10+I-1),I=1,6)
C # of backleg windings, # of turns and I for each winding   
      READ(NDAT,*) NBLWG,(A(NOEL,20+I),I=1,2*NBLWG)
      IF(NBLWG.GT.2) STOP ' SBR RAGSMM. # of blwg cannot exceed 2'
      A(NOEL,20) = NBLWG
C----- CHP FUITE ENTREE
C        XE, LambdaE, ff2,ff3
      READ(NDAT,*) (A(NOEL,I),I=30,33)
      READ(NDAT,*) II,(A(NOEL,40+I),I=1,II)
      A(NOEL,40) = II
C----- CHP FUITE SORTIE
      READ(NDAT,*) (A(NOEL,I),I=50,53)
      READ(NDAT,*) II,(A(NOEL,60+I),I=1,II)
      A(NOEL,60) = II
C----- Rotation of multipole components
      READ(NDAT,*) (A(NOEL,70+I-1),I=1,3)

C----- Step
      ND = 80
      CALL STPSIZ(NDAT,NOEL,ND,
     >                         A)

C KPOS=3 : XCE, YCE, ALE
C KPOS=4 : X-shft, Y-shft, Z-rot, Z-shft, Y-rot
      READ(NDAT,fmt='(a)') txt
      CALL STRGET(TXT,MSS,
     >                    NS,STRA) 
      
      READ(stra(1),*,err=66) II
      A(NOEL,90) = II 

      do kk = 2, NS
        READ(stra(kk),*,err=66) temp
        A(NOEL,91+kk-2) = temp
      enddo

      RETURN

 66   Continue
      RETURN
      END
