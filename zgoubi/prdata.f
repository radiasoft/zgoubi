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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE PRDATA(NLIN,FLIN,FDAT,
     >                                 LABEL,NOEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) FLIN, FDAT
      INCLUDE 'MXLD.H'
      CHARACTER(*) LABEL(MXL,*)
C-----------------------------------------------------------------------
C     Copy zgoubi.dat into zgoubi.res
C-----------------------------------------------------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG

      PARAMETER (I2000=2000)
      CHARACTER(I2000) TEXT
      PARAMETER (I6=6)
      CHARACTER(I6) TXT6
      INTEGER DEBSTR, FINSTR
      LOGICAL EMPTY
      LOGICAL IDLUNI, OK

      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LAB2(2)
      PARAMETER (KSIZ=10)
      PARAMETER (I104=104)
      DIMENSION LUNR(10)
      PARAMETER(MXFIL=1)
      CHARACTER(80) FINC(MXFIL)
      CHARACTER(50) CMMND

      WRITE(6,*) '  Copying  zgoubi.dat  into  zgoubi.res,'
      WRITE(6,*) '  numbering  and  labeling  elements...'

      OK = IDLUNI(
     >            NTMP)
      OPEN(UNIT=NTMP,FILE='zgoubi_temp.dat')

      IDA = 1
      LUNR(IDA) = NLIN
      LRD = LUNR(IDA)

C----- Read zgoubi.dat title (1st data line)
      READ(LRD,FMT='(A)',ERR=10,END=95) TEXT
      WRITE(NRES,FMT='(A)') TEXT(DEBSTR(TEXT):FINSTR(TEXT))
      WRITE(NTMP,FMT='(A)') TEXT(DEBSTR(TEXT):FINSTR(TEXT))

      NOEL=0
 10   CONTINUE
        READ (LRD,FMT='(A)',ERR=10,END=95) TEXT
        IF( .NOT. EMPTY(TEXT) ) THEN
          TEXT = TEXT(DEBSTR(TEXT):FINSTR(TEXT))
        ELSE
          TEXT = ' '
        ENDIF

        IDEB = 1
        IF( TEXT(IDEB:IDEB) .EQ. '''' ) THEN
          NOEL=NOEL+1
          TEXT = TEXT(DEBSTR(TEXT):I104)

          DO I=IDEB+1,IDEB+KSIZ+1
            IF(TEXT(I:I) .EQ. '''') GOTO 2
          ENDDO

 2        CONTINUE
          
          LABEL(NOEL,1) = ' '
          LABEL(NOEL,2) = ' '
          IF( .NOT. EMPTY(TEXT((I+1):I2000)) ) THEN
            CALL STRGET(TEXT((I+1):I2000),2,
     >                                      NST,LAB2)
            IF(NST.GE.1) THEN
              IF(LAB2(1)(1:1).NE.'!') THEN
                IF(.NOT. EMPTY(LAB2(1))) LABEL(NOEL,1) = LAB2(1)
                IF(NST.EQ.2) THEN
                  IF(LAB2(2)(1:1).NE.'!') THEN
                    IF(.NOT. EMPTY(LAB2(2))) LABEL(NOEL,2) = LAB2(2)
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF

          IF(TEXT(IDEB+1:I-1) .EQ. 'INCLUDE') THEN

            IF(FLIN(FINSTR(FLIN)-9:FINSTR(FLIN)) .EQ. 'zgoubi.dat')
     >      CALL ENDJOB('Pgm prdata. Job includes INCLUDE keyword, '//
     >      'hence name cannot be "zgoubi.dat" (reserved). Please give'
     >      //' different name. Use "zgoubi -fileIn filename" command.'
     >      ,-99)

            READ(LRD,FMT='(A)',ERR=10,END=95) TEXT
            READ(TEXT,*) NBFIL

            IF(NBFIL.GT.MXFIL) CALL ENDJOB('Pgm prdata. '//
     >      'INCLUDE has too many files. Max allowed is ',MXFIL)
            DO IFL = 1, NBFIL
              READ(LRD,FMT='(A)',ERR=10,END=95) FINC(IFL)
            ENDDO
            IDA = IDA + 1
            OK = IDLUNI(
     >                  LUNR(IDA)) 
            LRD = LUNR(IDA)
            OPEN(UNIT=LRD,FILE=FINC(1))           
            nOEL = NOEL -1 
            GOTO 10
          ELSE

            WRITE(TXT6,FMT='(I6)') NOEL
            TEXT = TEXT(1:I104)//TXT6
            WRITE(NRES,FMT='(T2,A)') TEXT(1:110)
            WRITE(NTMP,FMT='(T2,A)') TEXT(1:110)
          ENDIF

          IF(   TEXT(IDEB:IDEB+4) .EQ. '''FIN'''
     >     .OR. TEXT(IDEB:IDEB+4) .EQ. '''END''') GOTO 95

        ELSE

          WRITE(NRES,FMT='(A)') TEXT(DEBSTR(TEXT):FINSTR(TEXT))
          WRITE(NTMP,FMT='(A)') TEXT(DEBSTR(TEXT):FINSTR(TEXT))

        ENDIF

      GOTO 10

 95   CONTINUE
C      REWIND(LUNR(IDA))
      IF(IDA.GT.1) THEN
        CLOSE(LUNR(IDA))
        IDA = IDA-1
        LRD = LUNR(IDA)
        GOTO 10
      ENDIF

      CLOSE(NTMP)

      CMMND = 'mv zgoubi_temp.dat '//FDAT(DEBSTR(FDAT):FINSTR(FDAT))
      CALL SYSTEM(CMMND)

      RETURN
      END
