C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory 
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE PRDATA(
     >                  LABEL,NOEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'MXLD.H'
      CHARACTER LABEL(MXL,*)*8
C-----------------------------------------------------------------------
C     Copy zgoubi.dat into zgoubi.res
C-----------------------------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG

C      CHARACTER   TEXT*80
      PARAMETER (I110=110)
      CHARACTER   TEXT*110, txt6*6
      INTEGER DEBSTR, FINSTR
      LOGICAL EMPTY

      CHARACTER LAB2(2)*8
      PARAMETER (KSIZ=10)

      WRITE(6,*) '  Copying  zgoubi.dat  into  zgoubi.res,'
      WRITE(6,*) '  numbering  and  labeling  elements...'

C----- Read zgoubi.dat title (1st data line)
      READ (NDAT,FMT='(A)',ERR=10,END=95) TEXT
      IF(NRES .GT. 0) WRITE(NRES,FMT='(A)') TEXT

      NOEL=0
 10   CONTINUE
        READ (NDAT,FMT='(A)',ERR=10,END=95) TEXT
        IF( .NOT. EMPTY(TEXT) ) THEN
          TEXT = TEXT(DEBSTR(TEXT):FINSTR(TEXT))
        ELSE
          TEXT = ' '
        ENDIF
        IDEB = 1
        IF( TEXT(IDEB:IDEB) .EQ. '''' ) THEN
          NOEL=NOEL+1

          DO I=IDEB+1,IDEB+KSIZ+1
            IF(TEXT(I:I) .EQ. '''') GOTO 2
          ENDDO

 2        CONTINUE

            LABEL(NOEL,1) = '        '
            LABEL(NOEL,2) = '        '
            IF( .NOT. EMPTY(TEXT((I+1):I110)) ) THEN
              CALL STRGET(TEXT((I+1):I110),2,
     >                                       NST,LAB2)
              IF(LAB2(1)(1:1).NE.'!') THEN
                LABEL(NOEL,1) = LAB2(1)
                IF(NST.EQ.2) THEN 
                  IF(LAB2(2)(1:1).NE.'!') LABEL(NOEL,2) = LAB2(2)
                ENDIF
              ENDIF
            ENDIF

          WRITE(TXT6,FMT='(I6)') NOEL
          TEXT = TEXT(1:104)//TXT6
          IF(NRES .GT. 0) WRITE(NRES,FMT='(T2,A)') text

          IF(   TEXT(IDEB:IDEB+4) .EQ. '''FIN'''
     >     .OR. TEXT(IDEB:IDEB+4) .EQ. '''END''') GOTO 95

        ELSE
          IF(NRES .GT. 0) WRITE(NRES,FMT='(A)')  TEXT

        ENDIF

        GOTO 10

 95   CONTINUE
      REWIND(NDAT)

      RETURN
      END
