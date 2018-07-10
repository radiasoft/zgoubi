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
C  Brookhaven National Laboratory                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE RFIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ***************************************
C     READS DATA FOR FIT PROCEDURE WITH 'FIT'
C     ***************************************
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/MIMA/ DX(20)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/VARY/NV,IR(20),NC,I1(20),I2(20),V(20),IS(20),W(20),IC(20)
     > ,I3(20),XCOU(20)
 
      READ(NDAT,*) NV
      IF(NV.LT.1) RETURN
      DO 3 I=1,NV
 3      READ(NDAT,*) IR(I),IS(I),XCOU(I),DX(I)
      READ(NDAT,*) NC
      IF(NC.LT.1) RETURN
      DO 4 I=1,NC
 4      READ(NDAT,*) IC(I),I1(I),I2(I),I3(I),V(I),W(I)
      RETURN
      END
