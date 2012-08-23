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
      SUBROUTINE RBINAR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     *********************
C     READS DATA FOR BINARY
C     *********************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      PARAMETER (MXTA=40)
      COMMON/DONT/ TA(MXL,MXTA)

      CHARACTER*132 TXT132
      CHARACTER(80) STRA(2)

      DATA MXFIL / 9 /

C nf, ncol, NHdr : # of files to translate, # of data columns in file, # of header lines
      READ(NDAT,*) A(NOEL,1), A(NOEL,2), A(NOEL,3)
      NF = INT(A(NOEL,1))
      IF(NF.GT.20) CALL ENDJOB('Too many files, max is ',MXFIL)

C     ... FILE NAMES
      DO 1 I=1,NF
        READ(NDAT,200) TA(NOEL,I)
 200    FORMAT(A)
 1    CONTINUE
 
      RETURN
      END
