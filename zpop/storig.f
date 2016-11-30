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
      SUBROUTINE STORIG(NOEL,S,Y,Z,TETA1,XI,PHI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'MXLD.H'
      COMMON/LABCO/ ORIG(MXL,6)     
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) LABEL
      COMMON /LABEL/ LABEL(MXL,2)
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLEY
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ)
      INCLUDE 'MXFS.H'
      COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF),JPA(MXF,MXP)
      ORIG(NOEL,1)=S
      ORIG(NOEL,2)=Y
      ORIG(NOEL,3)=Z
      ORIG(NOEL,4)=TETA1
      ORIG(NOEL,5)=XI
      ORIG(NOEL,6)=PHI
      CALL FBGTXT
      WRITE(*,FMT='(A38,3A,1X,I6,1X,1P,5E14.6)') 
     >'KLEY,LBL1,LBL2,#Lmnt,S,Y,Z,TETA1,XI : '
     > ,KLEY,LABEL(NOEL,1),LABEL(NOEL,2),
     > NOEL,S,Y,Z,TETA1,XI
      RETURN
      END 
