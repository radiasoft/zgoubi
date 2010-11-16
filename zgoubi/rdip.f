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
C  Brookhaven National Laboratory               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE RDIP(ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     *********************
C     READS DATA FOR DIPOLE
C     *********************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
 
      NP = 1                                     ! IL      
      READ(NDAT,*) A(NOEL,NP)
      READ(NDAT,*) (A(NOEL,NP+I),I=1,2)          ! AT, RM
      NP=NP+2

C Magnet 1
      READ(NDAT,*) (A(NOEL,NP+I),I=1,5)          ! ACENT, HNORM, indices
      NP=NP+5

C       ... Entrance face
      READ(NDAT,*) (A(NOEL,NP+I),I=1,2)          ! LAMBDA, QSI
      NP=NP+2 
      READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! NBCOEF, COEFS_C0-5, SHIFT 
      NP=NP+8
      READ(NDAT,*) (A(NOEL,NP+I),I=1,6)          ! OMEGA,THETA,R1,U1,U2,R2
      NP=NP+6
C       ... Exit face 
      READ(NDAT,*) (A(NOEL,NP+I),I=1,2)          ! LAMBDA, QSI
      NP=NP+2 
      READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! NBCOEF, COEFS_C0-5, SHIFT 
      NP=NP+8
      READ(NDAT,*) (A(NOEL,NP+I),I=1,6)          ! OMEGA,THETA,R1,U1,U2,R2
      NP=NP+6
C       ... Lateral face
      READ(NDAT,*) (A(NOEL,NP+I),I=1,2)          ! LAMBDA, QSI
      NP=NP+2 
      READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! NBCOEF, COEFS_C0-5, SHIFT 
      NP=NP+8
      READ(NDAT,*) (A(NOEL,NP+I),I=1,7)          ! OMEGA,THETA,R1,U1,U2,R2,RM3
      NP=NP+7

C     ... IRD, RESOL
      READ(NDAT,*) A(NOEL,NP+1),A(NOEL,NP+2)
      NP=NP+2
C     ... XPAS
      NP=NP+1
      ND=NP
      READ(NDAT,*) A(NOEL,ND)
C     ... KP, RE, TE, RS, TS
C Modif, FM, Dec. 05
C      READ(NDAT,*) (A(NOEL,NP+I),I=1,5)
      READ(NDAT,*) (A(NOEL,NP+I),I=3,7)

      RETURN
      END
