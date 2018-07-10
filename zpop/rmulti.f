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
C  USA
C  -------
      SUBROUTINE RMULTI(NDAT,NOEL,MXL,A,MPOL,ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MXL,*)
C     ------------------------
C     READS DATA FOR MULTIPOLE
C     ------------------------
 
      CHARACTER(80) TXT*80

      IA = 1
      READ(NDAT,*) A(NOEL,IA) 

      IA = IA+1                
      IB = IA+MPOL+1           
      READ(NDAT,*) (A(NOEL,I),I=IA,IB)

      IA = IB+1                
C----- CHP FUITE ENTREE
C        XE, LambdaE, B1-Bmpol
      IB = IA + MPOL           
      READ(NDAT,*,ERR=11) (A(NOEL,I),I=IA,IB)
 11   CONTINUE
      IA = IB+1                
      IB = IA + 6              
      READ(NDAT,*) II,(A(NOEL,I),I=IA+1,IB)
      A(NOEL,IA) = II
C----- CHP FUITE SORTIE
      IA = IB + 1              
      IB = IA + MPOL           
      READ(NDAT,*,ERR=12) (A(NOEL,I),I=IA,IB)
 12   CONTINUE
      IA = IB+1                
      IB = IA + 6              
      READ(NDAT,*) II,(A(NOEL,I),I=IA+1,IB)
      A(NOEL,IA) = II
C----- Rotation of multipole components
      IA = IB +1               
      IB = IA + MPOL - 1       
      READ(NDAT,*) (A(NOEL,I),I=IA,IB)

      ND = IB + 1              
      READ(NDAT,FMT='(A)') TXT

      READ(NDAT,*) II,(A(NOEL,I),I=ND+4,ND+6)
      A(NOEL,ND+3) = II
 
 
      RETURN
      END
