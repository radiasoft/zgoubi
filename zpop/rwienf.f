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
      SUBROUTINE RWIENF(ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ************************
C     READS DATA FOR SEPARATOR
C     ************************
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      INCLUDE 'MXLD.H'        
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
 
C----- IL
      READ(NDAT,*) A(NOEL,1)
C----- XL(m), E(V/m), B(T), H-V
      READ(NDAT,*) (A(NOEL,I),I=10,13)
 
C----- ENTREE
C         XE,LMBD-E, LMBD-B
      READ(NDAT,*) (A(NOEL,I),I=20,22)
C         COEFS E
      READ(NDAT,*) (A(NOEL,I),I=30,35)
C         COEFS B
      READ(NDAT,*) (A(NOEL,I),I=40,45)
C----- SORTIE
C         XE,LMBD-E, LMBD-B
      READ(NDAT,*) (A(NOEL,I),I=50,52)
C         COEFS E
      READ(NDAT,*) (A(NOEL,I),I=60,65)
C         COEFS B
      READ(NDAT,*) (A(NOEL,I),I=70,75)
 
      ND = 80
C----- XPAS
      READ(NDAT,*) A(NOEL,80)
C----- KP,XCE,YCE,ALE
      READ(NDAT,*) A(NOEL,81),(A(NOEL,I),I=82,84)
 
      RETURN
      END
