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
      PROGRAM ZGOUBI_MAIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER * 9   DMY,HMS
      LOGICAL IDLUNI, READAT, FITING
      CHARACTER*10 FNAME

      PARAMETER (I5=5, I6=6)

      CALL INIDAT
      CALL RESET

      IF(IDLUNI(
     >          NDAT)) THEN
        FNAME = 'zgoubi.dat'
        OPEN(UNIT=NDAT,FILE=FNAME,STATUS='OLD',ERR=996)
      ELSE
        GOTO 996
      ENDIF

      IF(IDLUNI(
     >          NRES)) THEN
        FNAME = 'zgoubi.res'
        OPEN(UNIT=NRES,FILE=FNAME,ERR=997)
      ELSE
        GOTO 997
      ENDIF

      IF(IDLUNI(
     >          NLOG)) THEN
        FNAME = 'zgoubi.log'
        OPEN(UNIT=NLOG,FILE=FNAME,ERR=995)
      ELSE
        GOTO 995
      ENDIF

      CALL DATE2(DMY)
      CALL TIME2(HMS)
      CALL CPU_TIME(TIMSEC)

      WRITE(6   ,103) DMY,HMS
 103  FORMAT(/,'  Zgoubi, version 5.0.0.',/,
     >       '  Job  started  on  ',A,',  at  ',A)

      READAT = .TRUE.
      FITING = .FALSE.
      CALL FITSTA(I6,FITING)
      NL1 = 1
      NL2 = MXL
      CALL ZGOUBI(NL1,NL2,READAT)
      NOELMX=NL2

      CALL FITSTA(I5,
     >               FITING)
      IF(FITING) THEN
        READAT = .FALSE.
        CALL FITNU(*99)
        FITING = .FALSE.
        CALL FITSTA(I6,
     >                 FITING)
        NL2 = NOELMX-2   ! FIT keyword is at position NOELMX
C        write(*,*) ' zgoubi_main 2 fiting :',nres,nl1,nl2
        WRITE(6,201)
        WRITE(NRES,201)
        WRITE(6,200) 
        WRITE(NRES,200) 
 200    FORMAT(/,10X,
     >   ' MAIN PROGRAM :  now final run using FIT values ',A10)
        CALL ZGOUBI(NL1,NL2,READAT)
        WRITE(6,201)
        WRITE(NRES,201)
 201    FORMAT(/,128(1H*))
      ENDIF

      GOTO 10
 
 995  WRITE(6,*) ' PGM ZGOUBI : error open file ', FNAME
      GOTO 10
 996  WRITE(6,*) ' PGM ZGOUBI : error open file ', FNAME
      GOTO 10
 997  WRITE(6,*) ' PGM ZGOUBI : error open file ', FNAME
      GOTO 10 

 99     CONTINUE
        WRITE(NRES,103) DMY,HMS
        WRITE(6   ,103) DMY,HMS
        WRITE(NRES,101)
 101    FORMAT(/,128(1H*))
        WRITE(NRES,FMT='(/,10X,
     >          '' Main program : stopped upon key  FIT'')')

 10   CONTINUE

      WRITE(NRES,103) DMY,HMS
      WRITE(6   ,103) DMY,HMS
      CALL DATE2(DMY)
      CALL TIME2(HMS)
      WRITE(NRES,107) DMY,HMS
      WRITE(6   ,107) DMY,HMS
 107  FORMAT('  Job  ended  on  ',A,',  at  ',A,/)

      TEMP = TIMSEC
      CALL CPU_TIME(TIMSEC)
      WRITE(NRES,*) '  CPU time, total :  ',  TIMSEC-TEMP
      WRITE(   6,*) '  CPU time, total :  ',  TIMSEC-TEMP

      STOP
      END
