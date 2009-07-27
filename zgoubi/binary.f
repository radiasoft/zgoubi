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
C  Fran�ois M�ot <meot@lpsc.in2p3.fr>
C  Service Acc�lerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE BINARY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,20)
 
      CHARACTER*80 OLDFIL, NEWFIL
      LOGICAL BINARI, IDLUNI
      INTEGER DEBSTR,FINSTR
      DIMENSION X7(7)
      PARAMETER (I20=20)

      NFIC = INT(A(NOEL,1))
      NCOL = NINT(10.D0*A(NOEL,1)) - 10*NFIC

      IF(NCOL.LE.0) NCOL = 6

      IF(NFIC.GT.I20) 
     >   CALL ENDJOB('SBR  BINARY:  too  many  files,  max  is',I20)
 
      DO 1 IFIC=1,NFIC
        OLDFIL=TA(NOEL,IFIC)
C 200    FORMAT(A)
        IDEB = DEBSTR(OLDFIL)
        IFIN = FINSTR(OLDFIL)
 
        IF( BINARI(OLDFIL,IB) ) THEN
 
C          NEWFIL=OLDFIL(IDEB:IB-1)//OLDFIL(IB+2:IFIN)
          NEWFIL=OLDFIL(IB+2:IFIN)
          IF(IDLUNI(
     >              LNR)) THEN
            OPEN( UNIT=LNR, FILE=OLDFIL, FORM='UNFORMATTED'
     >      ,STATUS='OLD', ERR=96)
          ELSE
            GOTO 96
          ENDIF

          IF(IDLUNI(
     >              LNW)) THEN
C            OPEN( UNIT=LNW, FILE=NEWFIL, STATUS='NEW', ERR=97)
            OPEN( UNIT=LNW, FILE=NEWFIL, ERR=97)
          ELSE
            GOTO 97
          ENDIF

          IF(NRES.GT.0) THEN
            WRITE(NRES,100) OLDFIL, NEWFIL, NCOL
 100        FORMAT(10X,' Translate  from  binary  file  : ',A
     >          ,/,10X,' to  formatted  file            : ',A,/
     >          ,/,10X,' Number of data columns  : ',I4)
C            WRITE(6,200) ' TRANSLATING  FILE  ',OLDFIL
          ENDIF
 
          line = 0
 10       CONTINUE
            READ (LNR, ERR=90, END= 1 ) (X7(I),I=1,NCOL)
            line = line + 1
            WRITE(LNW,405) (X7(I),I=1,7)
 405        FORMAT(1X,1P,7E11.4)
          GOTO 10
 
        ELSE
 
C          NEWFIL=OLDFIL(IDEB:IB-1)//'B_'//OLDFIL(IB+2:IFIN)
          NEWFIL='B_'//OLDFIL(IB:IFIN)
          IF(IDLUNI(
     >              LNR)) THEN
            OPEN( UNIT=LNR, FILE=OLDFIL, STATUS='OLD', ERR=96)
          ELSE
            GOTO 96
          ENDIF
          IF(IDLUNI(
     >              LNW)) THEN
C            OPEN( UNIT=LNW, FILE=NEWFIL, FORM='UNFORMATTED'
C     >      , STATUS='NEW', ERR=97)
            OPEN( UNIT=LNW, FILE=NEWFIL, FORM='UNFORMATTED', ERR=97)
          ELSE
            GOTO 97
          ENDIF

 
          IF(NRES.GT.0) THEN
            WRITE(NRES,101) OLDFIL, NEWFIL, NCOL
 101        FORMAT(10X,' TRANSLATE  FROM  FORMATTED  FILE  : ',A
     >          ,/,10X,' TO  BINARY  FILE                  : ',A,/
     >          ,/,10X,' Number of data columns  : ',I4)
C            WRITE(6,200) ' TRANSLATING  FILE  ',OLDFIL
          ENDIF
 
          line = 0
 11       CONTINUE
C            READ (LNR,404, ERR=90, END= 1 ) (X7(I),I=1,NCOL)
            READ (LNR,*, ERR=90, END= 1 ) (X7(I),I=1,NCOL)
            line = line + 1
C 404        FORMAT(1X,7E11.2)
            WRITE(LNW) (X7(I),I=1,7)
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
      WRITE(NRES,*) ' Number of last line read is : ',line
      WRITE(6,*) ' Number of last line read is : ',line
      CALL ENDJOB('ERROR  DURING  READ  IN  '//OLDFIL,-99)
      RETURN
      END
