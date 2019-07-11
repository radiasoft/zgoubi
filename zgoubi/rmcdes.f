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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE RMCDES
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(132) TXT

      INTEGER DEBSTR

      READ(NDAT,FMT='(A)') TXT
      TXT = TXT(DEBSTR(TXT):LEN(TXT))
      IF(TXT(1:4).EQ.'INFO') THEN
        A(NOEL,4)=1.D0
        TXT = TXT(5:(LEN(TXT)-4))
      ELSE
        A(NOEL,4)=0.D0
      ENDIF

C----- Read M1, M2
      READ(TXT,*,ERR=99,END=99) (A(NOEL,I),I=1,2)
C----- Attempts reading M1, M2 and M1 life-time
      READ(TXT,*,ERR=10,END=10) (A(NOEL,I),I=1,3)
      GOTO 11
 10   CONTINUE
      A(NOEL,3) = 0.D0
C----- Seeds
 11   READ(NDAT,*) (A(NOEL,I),I=10,12)

      RETURN
 99   STOP ' *** DATA ERROR : in MCDESINT, while reading M1, M2 ***'
      END
