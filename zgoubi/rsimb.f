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
      SUBROUTINE RSIMB(ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************
C     READS DATA VENUS, PS170, QUADISEX, SEXQUAD
C     ******************************************
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/MARK/ KART,KALC,KERK,KUASEX
 
C     ... IL
      READ(NDAT,*) IA
      A(NOEL,1) = IA

C     ... XL, YL/RO, BO
      READ(NDAT,*) A(NOEL,10),A(NOEL,11),A(NOEL,12)

      ND=20

      IF(KUASEX .EQ. 3 .OR. KUASEX .EQ. 4.OR. KUASEX .EQ. 7) THEN
C       ... N, B1-2, G1-2 (KUASEX = 1,4) OR  R1,R2,R3,R4,BID (KUASEX = 7)
        READ(NDAT,*) (A(NOEL,I),I=20,24)
        ND = ND+10
      ENDIF

C     ... XPAS
      READ(NDAT,*) A(NOEL,ND)

C     ... KP,XCE, YCE, ALE
      ND1=ND+10
      READ(NDAT,*) IA,(A(NOEL,I),I=ND1+1,ND1+3)
      A(NOEL,ND1)=IA
 
      RETURN
      END
