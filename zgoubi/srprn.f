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
      SUBROUTINE SRPRN(KPR,LUN,IMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      CHARACTER*80 TA
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      PARAMETER (MXTA=45)
      COMMON/DONT/ TA(MXL,MXTA)
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM

      IF(KPR .GT. 1) THEN
        IF(IPASS.GT.1) THEN
          IF(KPR*(IPASS/KPR) .NE. IPASS) RETURN
        ENDIF
      ENDIF
 
C      IF(LUN.NE.NRES) THEN
C        IF(IPASS .EQ. 1) CALL OPEN2('SRPRN',
C     >                                      NSYN,TA(NOEL,1))
C      ENDIF

      CALL RAYSY2(IMAX,LUN)
      RETURN
      END
