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
      SUBROUTINE REPLAT(ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------------------
C     READS DATA FOR ELCYLDEF
C     ----------------------------------

      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL

C------ IL
      READ(NDAT,*) A(NOEL,1)
C------ Angle (rad), Layout radius (m), E field (V/m) at layout radius, index
      READ(NDAT,*) (A(NOEL,I),I=10,13)
C     ... XE, LE (cm, cm) and entrance fringe field coeffs
      READ(NDAT,*) (A(NOEL,I),I=20,21)
      READ(NDAT,*) (A(NOEL,I),I=30,35)
C     ... XS, LS and exit fringe field coeffs
      READ(NDAT,*) (A(NOEL,I),I=40,41)
      READ(NDAT,*) (A(NOEL,I),I=50,55)

      ND = 60
C     ... PAS
      READ(NDAT,*) A(NOEL,60)
C     ... KP
      READ(NDAT,*) KP
      A(NOEL,70) = KP
      IF( KP .EQ. 2 ) THEN
C       ... RE, TE, RS, TS
        READ(NDAT,*) (A(NOEL,I),I=80,83)
      ELSE
C       ... DP
        READ(NDAT,*) A(NOEL,80)
      ENDIF

      RETURN
      END
