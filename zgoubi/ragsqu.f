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
      SUBROUTINE RAGSQU(NDAT,NOEL,MXL,A,ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MXL,*)
C     -----------------------
C     READS DATA FOR AGS QUAD
C     -----------------------
      INCLUDE "C.MARK.H"     ! COMMON/MARK/ KART,KALC,KERK,KUASEX
      PARAMETER(MPOL=10)
 
C----- IL
      READ(NDAT,*) IA
      A(NOEL,1) = IA

C XL, R0, I_winding_1, I_winding_2, I_winding_3, dI1, dI2, dI3
      READ(NDAT,*) (A(NOEL,I),I=10,17)

C     ... Entrance fringe field
      READ(NDAT,*) (A(NOEL,I),I=20,21)
      READ(NDAT,*) IA,(A(NOEL,I),I=31,36)
      A(NOEL,30) = IA
C     ... Exit fringe field 
      READ(NDAT,*) (A(NOEL,I),I=40,41)
      READ(NDAT,*) IA,(A(NOEL,I),I=51,56)
      A(NOEL,50) = IA
C     ... Roll angle
      READ(NDAT,*) A(NOEL,60)
 
      ND=70
      CALL STPSIZ(NDAT,NOEL,ND,
     >                         A)

      READ(NDAT,*) IA,(A(NOEL,I),I=ND+10+1,ND+10+3)
      A(NOEL,ND+10) = IA

      RETURN
      END
