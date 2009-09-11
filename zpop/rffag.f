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
      SUBROUTINE RFFAG(ND,NDAT,
     >    NMAG,AT,RM,ACN,OP,XIE,OM,XIS,RE,TE,RS,TS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     --------------------------
C     READS DATA FOR SPIRAL FFAG
C     --------------------------
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
 
      READ(NDAT,*) A(NOEL,1)               ! IL      
      NP = 1                 
      READ(NDAT,*) (A(NOEL,NP+I),I=1,3)    ! NMAG, AT, R0
      NMAG = NINT(A(NOEL,NP+1))                                
      AT = A(NOEL,NP+2)
      RM = A(NOEL,NP+3)
      NP=NP+3

      DO 1 IMAG = 1, NMAG

        READ(NDAT,*) (A(NOEL,NP+I),I=1,4)          ! ACENT, DR0, HNORM, K
        ACN = A(NOEL,NP+1)
        NP=NP+4
C       ... Entrance face
        READ(NDAT,*) (A(NOEL,NP+I),I=1,2)          ! LAMBDA=g, gap's k 
        NP=NP+2                                    !    .eq.0/.ne.0 for constant/g_0(R0/r)^k
        READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! NBCOEF, COEFS_C0-5, SHIFT 
        NP=NP+8
        READ(NDAT,*) (A(NOEL,NP+I),I=1,6)          ! OMEGA,XI, 4*dummies (unused)
        OP = A(NOEL,NP+1)
        XIE = A(NOEL,NP+2)
        NP=NP+6
C         ... Exit face 
        READ(NDAT,*) (A(NOEL,NP+I),I=1,2)          ! LAMBDA=g, gap's k 
        NP=NP+2 
        READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! NBCOEF, COEFS_C0-5, SHIFT 
        NP=NP+8
        READ(NDAT,*) (A(NOEL,NP+I),I=1,6)          ! OMEGA,XI, 4*dummies (unused)
        OM = A(NOEL,NP+1)
        XIS = A(NOEL,NP+2)
        NP=NP+6
C         ... Lateral face
        READ(NDAT,*) (A(NOEL,NP+I),I=1,2)          ! LAMBDA=g, gap's k 
        NP=NP+2 
        READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! NBCOEF, COEFS_C0-5, SHIFT 
        NP=NP+8
        READ(NDAT,*) (A(NOEL,NP+I),I=1,6)          ! OMEGA,THETA,R1,U1,U2,R2
        NP=NP+6

 1    CONTINUE

C     ... IRD2(=0/1 for deriv. anal/interpol), RESOL -> mesh size= XPAS/RESOL
      READ(NDAT,*) A(NOEL,NP+1),A(NOEL,NP+2)
      NP=NP+2
C     ... XPAS
      NP=NP+1
      ND=NP
      READ(NDAT,*) A(NOEL,ND)
C     ... KP, RE, TE, RS, TS
      READ(NDAT,*) (A(NOEL,NP+I),I=3,7)
      RE = A(NOEL,NP+4)
      TE = A(NOEL,NP+5)
      RS = A(NOEL,NP+6)
      TS = A(NOEL,NP+7)
      WRITE(*,*) ' AT,RM,RE,TE,RS,TS : ',AT,RM,RE,TE,RS,TS

      RETURN
      END
