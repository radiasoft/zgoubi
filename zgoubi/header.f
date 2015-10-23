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
      SUBROUTINE HEADER(NL,NW,N,BINARY,
     >                                 *)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL BINARY
      INTEGER DEBSTR, FINSTR

      INCLUDE 'FILHDR.H'

      IF(NW.GT.0) WRITE(NW,FMT='(/,10X,''Now reading file header  ('',
     >I1,'' lines) : '',/)') N

      IF(.NOT.BINARY) THEN
        READ(NL,FMT='(A)',ERR=99,END=99) HDRF(1)
        IF(NW.GT.0) WRITE(NW,FMT='(A)') 
     >  HDRF(1)(DEBSTR(HDRF(1)):FINSTR(HDRF(1)))
        READ(NL,FMT='(A)',ERR=99,END=99) HDRF(2)
      ELSE
        READ(NL,ERR=99,END=89) HDRF(1)
        IF(NW.GT.0) WRITE(NW,FMT='(A)') 
     >  HDRF(1)(DEBSTR(HDRF(1)):FINSTR(HDRF(1)))
        READ(NL,ERR=99,END=89) HDRF(2)
      ENDIF
      IF(NW.GT.0) WRITE(NW,FMT='(A)') 
     >HDRF(2)(DEBSTR(HDRF(2)):FINSTR(HDRF(2)))

      IF(N.GT.MXHDF) CALL ENDJOB(
     >'Pgm header. Too many header lines in file. Sould be .le.',MXHDF)
      IF(.NOT.BINARY) THEN
        DO I=3, N
           READ(NL,FMT='(A)',ERR=99,END=89) HDRF(I)
           IF(NW.GT.0) WRITE(NW,FMT='(A)') 
     >     HDRF(I)(DEBSTR(HDRF(I)):FINSTR(HDRF(I)))
        ENDDO
      ELSE
        DO I=3, N
           READ(NL,          ERR=99,END=89) HDRF(1)
           IF(NW.GT.0) WRITE(NW,FMT='(A)') 
     >     HDRF(1)(DEBSTR(HDRF(1)):FINSTR(HDRF(1)))
        ENDDO
      ENDIF

      IF(NW.GT.0) WRITE(NW,FMT='(10X,''Header reading went on ok, it '',
     >''seems. Now proceeding... '',/)') 

      RETURN
 89   CONTINUE
      WRITE(6 ,*) 
     >'SBR HEADER : END of file reached while reading data file header'
      IF(NW.GT.0) WRITE(NW,*) 
     >'SBR HEADER : END of file reached while reading data file header'
      RETURN 1

 99   CONTINUE
      WRITE(6,*) 
     >'*** SBR HEADER : READ-error while reading data file header'
      WRITE(6,*) '        ... Empty file ?'
      IF(NW.GT.0) WRITE(NW,*) 
     >'*** SBR HEADER : READ-error while reading data file header'
      IF(NW.GT.0) WRITE(NW,*) '        ... Empty file ?'
      RETURN 1
      END
