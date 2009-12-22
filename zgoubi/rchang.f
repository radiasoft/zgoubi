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
      SUBROUTINE RCHANG
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER TXT*132, TXT1*1

      INTEGER DEBSTR
      PARAMETER (MSR=9,MSR2=2*MSR)
      CHARACTER QSHRO(MSR)*(2), SSHRO(MSR2)*(30)

      LOGICAL OLD
      DATA OLD / .TRUE. /

      READ(NDAT,FMT='(A)',ERR=99,END=99) TXT
      TXT = TXT(DEBSTR(TXT):LEN(TXT))
      TXT1 = TXT(1:1)
      IF( TXT1.EQ.'X' .OR.
     >    TXT1.EQ.'Y' .OR. 
     >    TXT1.EQ.'Z'     ) THEN
C New style, x-, y-, z-shift or  x-, y-, z-rotation in arbitrary order
        CALL STRGET(TXT,MSR2,
     >                       NSR2,SSHRO)
        NSR = NSR2/2
        DO I=1,NSR
          QSHRO(I) = SSHRO(2*I-1)(1:2)
          READ(SSHRO(2*I),*) A(NOEL,I)
        ENDDO
        OLD = .FALSE.
      ELSE
C old style, x- and y-shift followed by z-rotation
        NSR = 3
        READ(TXT,*,ERR=99,END=99) (A(NOEL,I),I=1,NSR)
        QSHRO(1) = 'XS'
        QSHRO(2) = 'YS'
        QSHRO(3) = 'ZR'
        OLD = .TRUE.
      ENDIF

      IF(OLD) QSHRO(4) = 'OL'

      CALL CHREF2(NSR,QSHRO)
      RETURN

 99   CONTINUE
      CALL ENDJOB(' End of job while reading from ''CHANGREF''',-99)
      RETURN

      END
