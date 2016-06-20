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
      SUBROUTINE RMATRX(
C     >                  IORD,IFOC,KWR,KCPL)
     >                  IORD,IFOC,KWR,SCPL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) SCPL
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL

      CHARACTER(132) TXT
      INTEGER DEBSTR
      LOGICAL STRCON

C      CALL MATRIC(NINT(
C     >  A(NOEL,1)),NINT(A(NOEL,2)),NINT(A(NOEL,3)),NINT(A(NOEL,4)))

      LINE = 1
      READ(NDAT,FMT='(A)',ERR=90,END=90) TXT
      TXT = TXT(DEBSTR(TXT):LEN(TXT))
C----- Read IORD, IFOC
      READ(TXT,*,ERR=90,END=90) IORD, IFOC
      A(NOEL,1) = IORD
      A(NOEL,2) = IFOC

      IF    (STRCON(TXT,'!',
     >                      IS)) TXT = TXT(1:IS-1)
      IF    (STRCON(TXT,'PRINT',
     >                          IS)) THEN
        KWR = 1
      ELSE
        KWR = 0
      ENDIF
      A(NOEL,3) = KWR
      
      IF    (STRCON(TXT,'coupled',
     >                            IS)) THEN
        SCPL = 'coupled'
C        KCPL = 1
      ELSE
        SCPL = 'uncoupled'
C        KCPL = 0
      ENDIF
      A(NOEL,4) = KCPL

      RETURN

 90   CALL ENDJOB('*** Pgm rmatrx, keyword MATRIX : '//
     >'input data error.',-99)
      RETURN
      END
