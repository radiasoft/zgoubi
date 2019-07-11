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
      SUBROUTINE RGOTO(NOEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
C      PARAMETER (LNTA=132) ; CHARACTER(LNTA) TA
C      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)

      CHARACTER(LNTA) TXT
      INTEGER DEBSTR
      LOGICAL STRCON

      LINE = 1
      READ(NDAT,FMT='(A)',ERR=97,END=97) TXT
      TXT = TXT(DEBSTR(TXT):LEN(TXT))
      IF(STRCON(TXT,'!',
     >                  IS)) TXT = TXT(1:IS-1)
      TA(NOEL,1) = TXT

      LINE = 2
      READ(NDAT,FMT='(A)',ERR=97,END=97) TXT
      TXT = TXT(DEBSTR(TXT):LEN(TXT))
      IF(STRCON(TXT,'!',
     >                  IS)) TXT = TXT(1:IS-1)
      TA(NOEL,2) = TXT

      RETURN

 97   WRITE(NRES,*)
     >'Data error met while reading condition from data list.'
      WRITE(6,*)
     >'Data error met while reading condition from data list.'
      GOTO 90

 90   CALL ENDJOB('*** Pgm rgoto. Input data error, at line ',line)
      RETURN
      END
