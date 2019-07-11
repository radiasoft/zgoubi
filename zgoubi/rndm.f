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
      FUNCTION RNDM()
c--------------------------------------------------------------------------------
C-- For ifort compiler. Comment otherwise
C BEWARE BEWARE BEWARE BEWARE BEWARE BEWARE BEWARE BEWARE BEWARE BEWARE BEWARE
C If you compile w ifort and This is commented, that may cause zgoubi not to work
#ifdef __INTEL_COMPILER
         use ifport
#endif
C--------------------------------------------------------------------------------
C While not required, the recommend usage of the Intel Fortran portability
C library functions is to access these either by inserting a USE IFPORT in
C the calling program, or by including iflport.f90 from the INCLUDE directory
C of your compiler distribution as part of your program.
C Use *without* the USE IFPORT (or iflport.f90) in conjunction with high-level
C optimizations may cause programs to suffer a run-time segmentation fault
C associated with calling the portability function. This is not a compiler defect.
C-----------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL R
      R = RAND(0)     ! RAND(1) RESTART THE SERIES
      RNDM = DBLE(R)
      RETURN

      ENTRY RNDM2(IRI)
      ISI = IRI
      CALL SRAND(ISI)
      RNDM2 = DBLE(ISI)
      RETURN

      END


