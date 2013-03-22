C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory       
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE ROPTIO
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      COMMON/DONT/ TA(MXL,MXTA)
      LOGICAL STRCON
      CHARACTER(80) TXT80
      INTEGER DEBSTR, FINSTR
C NY = 0/1 = off/on. NBOP = # of options, NBOP lines should follow
      READ(NDAT,*) NY, NBOP
      A(NOEL,1) = NY
      A(NOEL,2) = NBOP
      IF(NBOP.GT.40)
     >CALL ENDJOB('SBR roptio : nmbr of options exceded ; max is ',40)
      DO I = 1, NBOP
        READ(NDAT,FMT='(A)') TXT80
        TXT80 = TXT80(DEBSTR(TXT80):FINSTR(TXT80))
        IF(STRCON(TXT80,'!',
     >                      IS)) TXT80 = TXT80(DEBSTR(TXT80):IS-1)
        TA(NOEL,I) = TXT80
      ENDDO
      RETURN
      END