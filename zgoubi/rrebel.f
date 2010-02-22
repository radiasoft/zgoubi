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
      SUBROUTINE RREBEL(LABEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ***************************************
C     READS DATA FOR FIT PROCEDURE WITH 'FIT'
C     ***************************************
      INCLUDE 'MXLD.H'
      CHARACTER LABEL(MXL,2)*8
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL

      CHARACTER TXT132*132, TXTA*8, TXTB*8
      INTEGER DEBSTR, FINSTR
      LOGICAL STRCON
      CHARACTER STRA(3)*30

      READ(NDAT,FMT='(A)') TXT132
      TXT132 = TXT132(DEBSTR(TXT132):FINSTR(132))
      I3 = 3
      CALL STRGET(TXT132,I3,
     >                      NSR,STRA)

      READ(STRA(1),*,ERR=98) A(NOEL,1)
      READ(STRA(2),*,ERR=98) A(NOEL,2)
      READ(STRA(3),*,ERR=98) A(NOEL,3)
      IA3 = NINT(A(NOEL,3))
      A(NOEL,3) = IA3

      IF(STRCON(STRA(3),'.',
     >                      II)) THEN
        READ(STRA(3)(II+1:FINSTR(STRA(3))),*,ERR=98,END=98) IOP
        IF(IOP .GE. 1) THEN
          IF(IOP .EQ. 1) THEN
C Get Label, deduce related NOEL : multi-pass tracking will loop over NOEL-REBELOTE   
            READ(TXT132,*) DUM,DUM,DUM,TXTA
            DO JJ = 1, NOEL
              IF(LABEL(JJ,1).EQ.TXTA) THEN
                A(NOEL,4) = JJ
                GOTO 12
              ENDIF
            ENDDO           
            A(NOEL,4) = 1
 12         CONTINUE
            A(NOEL,5) = NOEL
          ELSEIF(IOP .EQ. 2) THEN
C Get 2 Labels, deduce related NOELs : multi-pass tracking will loop over NOEL1-NOEL2
            READ(TXT132,*) DUM,DUM,DUM,TXTA,TXTB
            DO JJ = 1, NOEL
              IF(LABEL(JJ,1).EQ.TXTA) THEN
                A(NOEL,4) = JJ
                GOTO 11
              ENDIF
            ENDDO           
            A(NOEL,4) = 1
 11         CONTINUE
            DO JJ = 1, NOEL
              IF(LABEL(JJ,1).EQ.TXTB) THEN
                A(NOEL,5) = JJ
                GOTO 10
              ENDIF
            ENDDO           
            A(NOEL,5) = NOEL
          ELSE
            CALL ENDJOB(' SBR RREBEL, no such option IOP=',IOP)
          ENDIF
        ELSE
          A(NOEL,4) = 1
          A(NOEL,5) = NOEL
        ENDIF
      ELSE
        A(NOEL,4) = 1
        A(NOEL,5) = NOEL
      ENDIF

 10   CONTINUE
      RETURN

 98   CONTINUE
      CALL ENDJOB(' SBR RREBEL, wrong input data / element #',-99)
      RETURN

      END
