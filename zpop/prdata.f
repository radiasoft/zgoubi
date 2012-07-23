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
      SUBROUTINE PRDATA(LABEL,NOEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'MXLD.H'
      CHARACTER LABEL(MXL,2)*(*)
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN

      CHARACTER   TITRE*80 , TXT*90, TXT6*6
      INTEGER DEBSTR
      LOGICAL EMPTY

      PARAMETER (LBLSIZ=10)
      CHARACTER LAB2(2)*(LBLSIZ)

      WRITE(6,*) '  Wait... reading .dat file'

C----- Read zgoubi.dat title (1st data line)
      READ (NDAT,FMT='(A)',ERR=10,END=95) TITRE

      NOEL=0
 10   CONTINUE
        READ (NDAT,FMT='(A)',ERR=10,END=95) TITRE
        IDEB=DEBSTR(TITRE)
        IF( TITRE(IDEB:IDEB) .EQ. '''' ) THEN
          NOEL=NOEL+1

          DO 1 I=IDEB+1,IDEB+9
            IF(TITRE(I:I) .EQ. '''') GOTO 2
 1        CONTINUE

 2        CONTINUE

            IF( .NOT. EMPTY(TITRE((I+1):80)) ) THEN
              CALL STRGET(TITRE((I+1):80),2,
     >                                             NST,LAB2)
              LABEL(NOEL,1) = LAB2(1)
              IF(NST.EQ.2) THEN 
                LABEL(NOEL,2) = LAB2(2)
              ELSE
                LABEL(NOEL,2) = '        '
               ENDIF
            ELSE       
              LABEL(NOEL,1) = '        '
              LABEL(NOEL,2) = '        '
            ENDIF

          WRITE(TXT6,FMT='(I6)') NOEL
          TXT = TITRE//'   '//TXT6
        ELSE
          TXT = TITRE
        ENDIF

        IF(   TITRE(IDEB:IDEB+4) .EQ. '''FIN'''
     >   .OR. TITRE(IDEB:IDEB+4) .EQ. '''END''') GOTO 95
        GOTO 10

 95   CONTINUE
      REWIND(NDAT)
      RETURN
      END
