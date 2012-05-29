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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE REVERS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     Copy zgoubi.dat into zgoubi.res
C-----------------------------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG

      CHARACTER   TEXT*110
      INTEGER DEBSTR
      LOGICAL IDLUNI
      CHARACTER*12 NAMSTO, NAMTMP
      PARAMETER (NAMSTO='reverse.dat',NAMTMP='reverse2.dat')

      WRITE(6,*) '  Reversing  zgoubi.dat  into  reverse.dat,'

          IF(IDLUNI(
     >              LUN)) THEN
              OPEN(UNIT=LUN,FILE=NAMSTO,ERR=99)
              CLOSE(UNIT=LUN,STATUS='DELETE')
              OPEN(UNIT=LUN,FILE=NAMSTO,ERR=99)
          ENDIF
          IF(IDLUNI(
     >              LUN2)) THEN
              OPEN(UNIT=LUN2,FILE=NAMTMP,ERR=99)
              CLOSE(UNIT=LUN2,STATUS='DELETE')
              OPEN(UNIT=LUN2,FILE=NAMTMP,ERR=99)
          ENDIF

      REWIND(NDAT)

C----- Read zgoubi.dat title (1st data line)
      LSTO = LUN
      LTMP = LUN2
      READ (NDAT,FMT='(A)',ERR=10,END=95) TEXT
      WRITE(LSTO,FMT='(A)') TEXT

      NOEL=0
 10   CONTINUE

        READ (NDAT,FMT='(A)',ERR=10,END=95) TEXT
        IDEB=DEBSTR(TEXT)

        IF( TEXT(IDEB:IDEB) .EQ. '''' ) THEN
          CALL APPEN(LTMP,LSTO)

              CLOSE(UNIT=LTMP,STATUS='DELETE')
              OPEN(UNIT=LTMP,FILE=NAMTMP,ERR=99)
          CALL APPEN(LSTO,LTMP)
              CLOSE(UNIT=LSTO,STATUS='DELETE')
              OPEN(UNIT=LSTO,FILE=NAMSTO,ERR=99)
        
          NOEL=NOEL+1

        ENDIF

          WRITE(LSTO,FMT='(A,I6)') TEXT, NOEL

          IF(   TEXT(IDEB:IDEB+4) .EQ. '''FIN'''
     >     .OR. TEXT(IDEB:IDEB+4) .EQ. '''END''') GOTO 95


        GOTO 10

 95   CONTINUE

          CALL APPEN(LTMP,LSTO)
              CLOSE(UNIT=LTMP,STATUS='DELETE')
        
      RETURN

 99   IF(NRES .GT. 0) WRITE(NRES,*) ' sbr revers, error open '
      RETURN
      END
