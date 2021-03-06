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
      SUBROUTINE RCYCLO(ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     --------------------------
C     READS DATA FOR SPIRAL CYCLOTRON
C     --------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL

      CHARACTER(132) TXT

      READ(NDAT,*) A(NOEL,1)               ! IL
      NP = 1
      READ(NDAT,*) (A(NOEL,NP+I),I=1,4)    ! NMAG, AT, R0, Type of sector (radial, spiral, both)
      NMAG = NINT(A(NOEL,NP+1))
      NP=NP+4



      DO 1 IMAG = 1, NMAG

        READ(NDAT,*) (A(NOEL,NP+I),I=1,10)          ! ACENT, DR0, FAC, HNORM, K, Rref, H1,H2,H3,H4
        NP=NP+10
C       ... Entrance face
        READ(NDAT,*) (A(NOEL,NP+I),I=1,4)          ! LAMBDA=g, gap's k  g10 g11
        NP=NP+4                                    !    .eq.0/.ne.0 for constant/g_0(R0/r)^k
        READ(NDAT,*) (A(NOEL,NP+I),I=1,10)          ! NBCOEF, COEFS_C0-7, NORME
        NP=NP+10

        READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! OMEGA,XI0,XI1,XI2,XI3,aen,ben,cen
        NP=NP+8
C         ... Exit face
        READ(NDAT,*) (A(NOEL,NP+I),I=1,4)          ! LAMBDA=g, gap's k g20 g21
        NP=NP+4
        READ(NDAT,*) (A(NOEL,NP+I),I=1,10)          ! NBCOEF, COEFS_C0-7, NORMS
        NP=NP+10
        READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! OMEGA,XI0exit,XI1exit,XI2exit,XI3exit,aexit,bexit,cexit
        NP=NP+8
C         ... Lateral face
        READ(NDAT,*) (A(NOEL,NP+I),I=1,2)          ! LAMBDA=g, gap's k
        NP=NP+2
        READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! NBCOEF, COEFS_C0-5, SHIFT
        NP=NP+8
        READ(NDAT,*) (A(NOEL,NP+I),I=1,6)          ! OMEGA,THETA,R1,U1,U2,R2
        NP=NP+6

 1    CONTINUE

C KIRD= 0 or 2,4,25  for analytic or 2-,4-,5-type numerical interpolation
C mesh size= XPAS/RESOL
C      READ(NDAT,*) A(NOEL,NP+1),A(NOEL,NP+2)
      READ(NDAT,FMT='(A132)') TXT
      READ(TXT,*,ERR=999) IA
      IF(IA.NE.0) THEN
C KIRD=2,25 OR 4 AND RESOL
        READ(TXT,*,ERR=999) IA,A(NOEL,NP+2)
        A(NOEL,NP+1) = IA
      ELSE
C KIRD=2,25 OR 4 AND IRD==RESOL
        READ(TXT,*,ERR=999) IA,IB
        A(NOEL,NP+1) = IA
        A(NOEL,NP+2) = IB
      ENDIF
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

 999  CONTINUE
      CALL ENDJOB('SBR rffag. Input data format error at KIRD.',-99)
      RETURN
      END
