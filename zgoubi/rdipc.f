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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE RDIPC(ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     **********************
C     READS DATA FOR DIPOLES
C     **********************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
 
      NP = 1                                     ! IL      
      READ(NDAT,*) A(NOEL,NP)
      READ(NDAT,*) (A(NOEL,NP+I),I=1,3)          ! NMAG, XTNT, RM
      NMAG = NINT(A(NOEL,NP+1))
      NP=NP+3

      DO 1 IMAG = 1, NMAG

        READ(NDAT,*) 
     >    (A(NOEL,NP+I),I=1,3),IND,(A(NOEL,NP+4+I),I=1,ABS(IND)) ! XCENT, DYM, HNORM, indices
        A(NOEL,NP+4) = IND
        NP=NP+4 + ABS(IND)

C       ... Entrance face
        READ(NDAT,*) (A(NOEL,NP+I),I=1,2)          ! LAMBDA, QSI
        NP=NP+2 
        READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! NBCOEF, COEFS_C0-5, SHIFT 
        NP=NP+8
        READ(NDAT,*) (A(NOEL,NP+I),I=1,6)          ! DXP,THETA,R1,U1,U2,R2
        NP=NP+6
C         ... Exit face 
        READ(NDAT,*) (A(NOEL,NP+I),I=1,2)          ! LAMBDA, QSI
        NP=NP+2 
        READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! NBCOEF, COEFS_C0-5, SHIFT 
        NP=NP+8
        READ(NDAT,*) (A(NOEL,NP+I),I=1,6)          ! DXM,THETA,R1,U1,U2,R2
        NP=NP+6
C         ... Lateral face
        READ(NDAT,*) (A(NOEL,NP+I),I=1,2)          ! LAMBDA, QSI
        NP=NP+2 
        READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! NBCOEF, COEFS_C0-5, SHIFT 
        NP=NP+8
        READ(NDAT,*) (A(NOEL,NP+I),I=1,7)          ! DXL,THETA,R1,U1,U2,R2,RM3
        NP=NP+7

 1    CONTINUE

C     ... IRD, RESOL
      READ(NDAT,*) A(NOEL,NP+1),A(NOEL,NP+2)
      NP=NP+2
C     ... XPAS
      NP=NP+1
      ND=NP
      READ(NDAT,*) A(NOEL,ND)
C     ... KP, XCE, YCE, ALE
      READ(NDAT,*) (A(NOEL,NP+I),I=3,6)

      RETURN
      END
