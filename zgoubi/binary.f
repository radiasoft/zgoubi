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
      SUBROUTINE BINARY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
 
      CHARACTER(80) OLDFIL, NEWFIL
      CHARACTER(132) HEADER
      LOGICAL BINARI, IDLUNI
      INTEGER DEBSTR,FINSTR
      DIMENSION X7(7)
      PARAMETER (I20=20)
      CHARACTER(40) FRMT

      NFIC = INT(A(NOEL,1))
      NFRM =NINT( 10*(A(NOEL,1) - INT(A(NOEL,1))))
      NCOL = NINT(A(NOEL,2))
      NHEAD = NINT(A(NOEL,3))
      IF(NFIC.GT.I20) 
     >   CALL ENDJOB('SBR  BINARY:  too  many  files,  max  is',I20)
      IF(NCOL.LE.0) NCOL = 6
C May be second argument in BINARY (zgoubi version > 5.1.0) : # of header lines
 
      DO 1 IFIC=1,NFIC
        OLDFIL=TA(NOEL,IFIC)
        OLDFIL=OLDFIL(DEBSTR(OLDFIL):FINSTR(OLDFIL))
        IDEB = 1
        IFIN = FINSTR(OLDFIL)
 
        IF( BINARI(OLDFIL,IB) ) THEN
 
          NEWFIL=OLDFIL(IB+2:IFIN)
          IF(NRES.GT.0) THEN
            WRITE(NRES,100) OLDFIL, NEWFIL, NCOL, NHEAD, NFRM
 100        FORMAT(10X,' Translate  from  binary  file  : ',A
     >          ,/,10X,' to  formatted  file            : ',A,/
     >          ,/,10X,' Number of data columns  : ',I4
     >          ,/,10X,' Number of header lines  : ',I4
     >          ,/,10X,' Format type  : ',I4)
          ENDIF
 
          IF(IDLUNI(
     >              LNR)) THEN
            OPEN( UNIT=LNR, FILE=OLDFIL, FORM='UNFORMATTED'
     >      ,STATUS='OLD', ERR=96)
          ELSE
            GOTO 96
          ENDIF

          IF(IDLUNI(
     >              LNW)) THEN
            OPEN( UNIT=LNW, FILE=NEWFIL, ERR=97)
          ELSE
            GOTO 97
          ENDIF

          LINE = 0
          DO NH = 1, NHEAD
            READ (LNR, ERR=90, END= 1 ) HEADER
            LINE = LINE + 1
            WRITE(LNW,FMT='(A)') HEADER(DEBSTR(HEADER):FINSTR(HEADER))
          ENDDO

 10       CONTINUE
            READ (LNR, ERR=90, END= 1 ) (X7(I),I=1,NCOL)
            LINE = LINE + 1
            WRITE(FRMT,FMT='(A,I0,A)') '(1P,',NCOL,'(1X,E20.12))'
            WRITE(LNW,FRMT) (X7(I),I=1,NCOL)

          GOTO 10
 
        ELSE
 
          NEWFIL='b_'//OLDFIL(IB:IFIN)
          IF(NRES.GT.0) THEN
            WRITE(NRES,101) OLDFIL, NEWFIL, NCOL, NHEAD, NFRM
 101        FORMAT(10X,' TRANSLATE  FROM  FORMATTED  FILE  : ',A
     >          ,/,10X,' TO  BINARY  FILE                  : ',A,/
     >          ,/,10X,' Number of data columns  : ',I4
     >          ,/,10X,' Number of header lines  : ',I4
     >          ,/,10X,' Format type  : ',I4)
C            WRITE(6,200) ' TRANSLATING  FILE  ',OLDFIL
          ENDIF
 
          IF(IDLUNI(
     >              LNR)) THEN
c               write(*,*) ' sbr binary, try open ',oldfil,lnr
            OPEN( UNIT=LNR, FILE=OLDFIL, ERR=96)
          ELSE
            WRITE(ABS(NRES),*) ' SBR BINARY, no idle unit for'
     >            ,' opening file ',oldfil
            GOTO 96
          ENDIF
          IF(IDLUNI(
     >              LNW)) THEN
            OPEN( UNIT=LNW, FILE=NEWFIL, FORM='UNFORMATTED', ERR=97)
          ELSE
            GOTO 97
          ENDIF
 
          LINE = 0
          DO NH = 1, NHEAD
            READ (LNR,FMT='(A132)', ERR=90, END= 1 ) HEADER
            LINE = LINE + 1
            WRITE(LNW) HEADER
          ENDDO
 11       CONTINUE
            IF    (NFRM.EQ.0) THEN
              READ (LNR,*, ERR=90, END= 1 ) (X7(I),I=1,NCOL)
            ELSEIF(NFRM.EQ.1) THEN
              WRITE(FRMT,FMT='(A,I0,A)') '(1X,',NCOL,'E11.2)'
              READ (LNR,FMT=FRMT) (X7(I),I=1,NCOL)
            ENDIF
            LINE = LINE + 1
            WRITE(LNW) (X7(I),I=1,NCOL)
          GOTO 11
        ENDIF
 
 1    CONTINUE
      CLOSE(LNR)
      CLOSE(LNW)
 
      IF(NRES .GT. 0) WRITE(NRES,102) NFIC
 102  FORMAT(//,15X,' TRANSLATION  OF ',I2,'  FILES  COMPLETED')
      RETURN
 
 96   CALL ENDJOB('ERROR  OPEN  OLD  FILE '//OLDFIL,-99)
 97   CALL ENDJOB('ERROR  OPEN  NEW  FILE '//NEWFIL,-99)
 90   CONTINUE
      WRITE(NRES,*) ' '
      WRITE(NRES,*) ' Number of last line read is : ',line
      CALL ENDJOB('ERROR  DURING  READ  IN  '//OLDFIL,-99)
      RETURN
      END
