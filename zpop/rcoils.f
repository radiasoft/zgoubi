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
C  Brookhaven National Laboratory                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE RCOILS(NDAT,NOEL,MXL,A,ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MXL,*)
C     ------------------------
C     READS DATA FOR COILS
C     ------------------------
      COMMON/MARK/ KART,KALC,KERK,KUASEX
 
C------ IL
      READ(NDAT,*) A(NOEL,1)
C MS= total number of solenoids, D2, ... , DMS = distances from center to center. 
C----- MS, XL_1, RO_1, BO_1, Dist_2, XL_2, RO_2, BO_2, ... , Dist_MS, XL_MS, RO_MS, BO_MS
      READ(NDAT,*) MS,
     >(A(NOEL,3+4*M),A(NOEL,4+4*M),A(NOEL,5+4*M),A(NOEL,6+4*M),M=0,MS-1)
      A(NOEL,2) = MS
C----- XE, XS
      IA = 7+4*MS
      READ(NDAT,*) A(NOEL,IA),A(NOEL,IA+1)

C----- Integr. step
      ND = IA+2

C      CALL STPSIZ(NDAT,NOEL,
C     >                      A,ND)

      READ(NDAT,*) JA,(A(NOEL,I),I=ND+2,ND+4)
      A(NOEL,ND+1) = JA
 
      RETURN
      END
