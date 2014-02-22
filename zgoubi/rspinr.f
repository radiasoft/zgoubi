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
      SUBROUTINE RSPINR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     **************************
C     READ DATA FOR SPIN ROTATOR
C     **************************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(132) TXT
      INTEGER DEBSTR
      LOGICAL STRCON

C On/Off switch
      READ(NDAT,*) IOP
      A(NOEL,1) = IOP

C Old form : ang=rotationAng.axisAng
C      READ(NDAT,*) ANG
C      A(NOEL,10) = ANG

      READ(NDAT,FMT='(A)') TXT
      IF(STRCON(TXT,'!',
     >                  IS)) TXT = TXT(DEBSTR(TXT):IS-1)
      READ(TXT,*,END=10,ERR=10) ANG, PHI
      A(NOEL,10) = ANG
      A(NOEL,11) = PHI
      RETURN

 10   WRITE(NRES,*) ' Format error in SPINR input data'
      CALL ENDJOB('SBR RSPINR. Data must have the form angle, phi',-99)
      RETURN
      END
