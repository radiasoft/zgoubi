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
      SUBROUTINE RSRLOS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,40)
      CHARACTER(80) STRA(2)

      READ(NDAT,*,ERR=99) A(NOEL,1)
      READ(NDAT,FMT='(A80)',ERR=99) TA(NOEL,1)

      CALL STRGET(TA(NOEL,1),2,
     >                         NOUT,STRA)

      TA(NOEL,1)=' '
      TA(NOEL,2)=' '
      IF(NOUT.GE.1) TA(NOEL,1)=STRA(1)
      IF(NOUT.GE.2) TA(NOEL,2)=STRA(2)

      READ(NDAT,*,ERR=99) A(NOEL,10), A(NOEL,11)

      RETURN
 99   CONTINUE
      WRITE(6,*) ' Input data error at keyword SRLOSS'
      STOP
      END
