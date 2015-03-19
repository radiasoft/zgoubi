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
      SUBROUTINE RFFAG(
     >                 ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     --------------------------
C     READS DATA FOR SPIRAL FFAG
C     --------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
 
      CHARACTER(400) TXT
      INTEGER DEBSTR

      READ(NDAT,*) A(NOEL,1)               ! IL      
      NP = 1                 
      READ(NDAT,*) (A(NOEL,NP+I),I=1,3)    ! NMAG, AT, R0
      NMAG = NINT(A(NOEL,NP+1))                                
      NP=NP+3

      DO IMAG = 1, NMAG

        READ(NDAT,*) (A(NOEL,NP+I),I=1,4)          ! ACENT, DR0, HNORM, K
        NP=NP+4
C       ... Entrance face
        READ(NDAT,*) (A(NOEL,NP+I),I=1,2)          ! LAMBDA=g, gap's k 
        NP=NP+2                                    !    .eq.0/.ne.0 for constant/g_0(R0/r)^k
        READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! NBCOEF, COEFS_C0-5, SHIFT 
        NP=NP+8
        READ(NDAT,*) (A(NOEL,NP+I),I=1,6)          ! OMEGA,XI, 4*dummies (unused)
        NP=NP+6
C         ... Exit face 
        READ(NDAT,*) (A(NOEL,NP+I),I=1,2)          ! LAMBDA=g, gap's k 
        NP=NP+2 
        READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! NBCOEF, COEFS_C0-5, SHIFT 
        NP=NP+8
        READ(NDAT,*) (A(NOEL,NP+I),I=1,6)          ! OMEGA,XI, 4*dummies (unused)
        NP=NP+6
C         ... Lateral face
        READ(NDAT,*) (A(NOEL,NP+I),I=1,2)          ! LAMBDA=g, gap's k 
        NP=NP+2 
        READ(NDAT,*) (A(NOEL,NP+I),I=1,8)          ! NBCOEF, COEFS_C0-5, SHIFT 
        NP=NP+8
        READ(NDAT,*) (A(NOEL,NP+I),I=1,6)          ! OMEGA,THETA,R1,U1,U2,R2
        NP=NP+6

      ENDDO

      READ(NDAT,FMT='(A400)') TXT
      TA(NOEL,1) = ' '

      IF(TXT(DEBSTR(TXT):DEBSTR(TXT)+5) .EQ. 'IntLim') THEN
C Integration limts defined by entrance and/or exit lines
        READ(TXT(DEBSTR(TXT)+6:400),*) IDRT  
        NP = NP + 1
        A(NOEL,NP) = IDRT
        IF    (IDRT.EQ.-1) THEN
C          Intgration limit at entrance
          READ(TXT(DEBSTR(TXT)+6:400),*) IDUM,(A(NOEL,NP+I),I=1,3)
          NP = NP + 3
        ELSEIF(IDRT.EQ. 1) THEN
C          Intgration limit at exit
          READ(TXT(DEBSTR(TXT)+6:400),*) IDUM,(A(NOEL,NP+I),I=1,3)
          NP = NP + 3
        ELSEIF(IDRT.EQ. 2) THEN
C          Intgration limit at entrance and at exit
          READ(TXT(DEBSTR(TXT)+6:400),*) IDUM,(A(NOEL,NP+I),I=1,6)
          NP = NP + 6
        ELSE
          CALL ENDJOB('Pgm rffag. No such option IDRT = ',IDRT)
        ENDIF
        TA(NOEL,1) = TXT(DEBSTR(TXT):DEBSTR(TXT)+5)
        READ(NDAT,FMT='(A400)') TXT
      ENDIF

C KIRD= 0 or 2,4,25  for analytic or 2-,4-,5-type numerical interpolation
C mesh size= XPAS/RESOL
C      READ(NDAT,*) A(NOEL,NP+1),A(NOEL,NP+2)
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
C      READ(NDAT,*) (A(NOEL,NP+I),I=3,7)
      READ(NDAT,*) (A(NOEL,NP+I),I=1,5)
      KP = NINT(A(NOEL,NP+1))

      RETURN

 999  CONTINUE
      CALL ENDJOB('SBR rffag. Input data format error at KIRD.',-99)
      RETURN
      END
