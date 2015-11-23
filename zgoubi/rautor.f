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
      SUBROUTINE RAUTOR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------
C     READS DATA FOR AUTREF
C     ---------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(132) TXT132
      LOGICAL STRCON
      INTEGER DEBSTR, FINSTR
      READ(NDAT,*) TXT132
      TXT132 = TXT132(DEBSTR(TXT132):FINSTR(TXT132))
      IF(STRCON(TXT132,'!',
     >                     ISA)) THEN
        TXT132 = TXT132(1:ISA-1)
      ELSE
        ISA = FINSTR(TXT132)
      ENDIF
      IF(STRCON(TXT132,'.',
     >                     IS)) THEN
        READ(TXT132(   1:IS-1),*) IA
        READ(TXT132(IS+1:ISA),*) IA2
      ELSE
        READ(TXT132,*) IA
        IA2 = 0
      ENDIF
      A(NOEL,1) = IA
      A(NOEL,2) = IA2
C Read 3 particle numbers
      IF    (IA .EQ. 3) THEN
        READ(NDAT,*) A(NOEL,10),A(NOEL,11),A(NOEL,12)
C Read 3 centering coordinates
      ELSEIF(IA .EQ. 4) THEN
        IF    (IA2.EQ.0) THEN
C Center the beam on xce, yce, ale
          READ(NDAT,*) A(NOEL,10),A(NOEL,11),A(NOEL,12)
        ELSEIF(IA2.LE.2) THEN
C 1 : Center the beam on xce, yce, ale, p/pRef, time
C 2 : Center the beam on xce, yce, ale, p/pRef and set all times to A(NOEL,14)
          READ(NDAT,*) A(NOEL,10),A(NOEL,11),A(NOEL,12),
     >    A(NOEL,13),A(NOEL,14)
        ELSE
          CALL ENDJOB('Pgm rautor. No such option IA2 = ',IA2)
        ENDIF
      ELSEIF(IA .EQ. 5) THEN
        READ(NDAT,*) A(NOEL,10),A(NOEL,11)
      ELSE
C        CALL ENDJOB('Pgm RAUTOR. No such option IA = ',IA2)
      ENDIF
      RETURN
      END
