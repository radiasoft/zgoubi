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
      PROGRAM ZGOUBI_MAIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(10) DMY
      CHARACTER(9) HMS
      LOGICAL IDLUNI, READAT, FITING, ENDFIT
      CHARACTER(100) FILIN, FILOU, FILOG
      CHARACTER(10) FDAT, FRES, FLOG

      PARAMETER (I5=5, I6=6)

      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLE
      LOGICAL FITFNL
      INTEGER DEBSTR, FINSTR

      character(len=12), dimension(:), allocatable :: args
      logical savxec, savzpp

      DATA ENDFIT / .FALSE. /
      DATA FDAT, FRES, FLOG / 'zgoubi.dat', 'zgoubi.res', 'zgoubi.log'/
      data savxec, savzpp / .false., .false.  /

C Manage possible arguments to zgoubi -----------------------
      NBARGS = COMMAND_ARGUMENT_COUNT()
      ALLOCATE(ARGS(NBARGS)) 
      DO IX = 1, NBARGS
         CALL GET_COMMAND_ARGUMENT(IX,ARGS(IX))
      END DO
      IF(NBARGS.GE.1) THEN
        WRITE(*,*) ' '
        WRITE(*,*) ' Argument(s) to zgoubi command :'
        WRITE(*,*) (IX,ARGS(IX),IX = 1, NBARGS)
      ENDIF
      DO IX = 1, NBARGS
         IF(ARGS(IX) .EQ. '-fileIn'  ) FDAT = ARGS(IX+1)
         IF(ARGS(IX) .EQ. '-fileOut' ) FRES = ARGS(IX+1)
         IF(ARGS(IX) .EQ. '-fileLog' ) FLOG = ARGS(IX+1)
         IF(ARGS(IX) .EQ. '-saveExec') SAVXEC = .TRUE.
         IF(ARGS(IX) .EQ. '-saveZpop') SAVZPP = .TRUE.
      ENDDO
C -----

      IF(IDLUNI(
     >          NDAT)) THEN
        FILIN = FDAT
        OPEN(UNIT=NDAT,FILE=FILIN,STATUS='OLD',ERR=996)
c          write(*,*) ' zgoubi.dat is unit # ',ndat
      ELSE
        GOTO 996
      ENDIF

      IF(IDLUNI(
     >          NRES)) THEN
        FILOU = FRES
        OPEN(UNIT=NRES,FILE=FILOU,ERR=997)
c          write(*,*) ' zgoubi.res is unit # ',nres
      ELSE
        GOTO 997
      ENDIF

      IF(IDLUNI(
     >          NLOG)) THEN
        FILOG = FLOG
        OPEN(UNIT=NLOG,FILE=FILOG,ERR=995)
      ELSE
        GOTO 995
      ENDIF

      IF(SAVXEC .OR. SAVZPP) THEN
        WRITE(6,*) ' '
        WRITE(6,*) ' --'
        if(savxec) then
          WRITE(6,fmt='(a)') 
     >    'Now copying zgoubi executable to local directory.'
          CALL SYSTEM('cp ~/zgoubi/SVN/current/zgoubi/zgoubi .') 
        endif
        if(savzpp) then
          WRITE(6,fmt='(a)') 'Now copying zpop executable to local '
     >    //'directory (if you wish to run it, use an xterm terminal).'
          CALL SYSTEM('cp ~/zgoubi/SVN/current/zpop/zpop .') 
        endif
        WRITE(6,*) ' --'
        WRITE(6,*) ' '
      ENDIF

      CALL DATE2(DMY)
      CALL TIME2(HMS)
      CALL CPU_TIME(TIMSEC)

      WRITE(6   ,103) DMY,HMS
 103  FORMAT(/,'  Zgoubi, author''s dvlpmnt version.',/,
     >       '  Job  started  on  ',A,',  at  ',A)

 1    CONTINUE

      READAT = .TRUE.

 11   CONTINUE
      FITING = .FALSE.
      CALL INIDAT
      CALL RESET
      CALL CHECKS
      CALL FITSTA(I6,FITING)
      NL1 = 1
      NL2 = MXL
      ENDFIT = .FALSE.
      CALL ZGOUBI(NL1,NL2,READAT,
     >                           NBLMN,ENDFIT)

      CALL FITSTA(I5,
     >               FITING)
C          write(*,*) ' main ',fiting,readat,endfit
C              pause
      IF(FITING) THEN
        READAT = .FALSE.
        CALL FITNU(NRES,*99)
        FITING = .FALSE.
        CALL FITSTA(I6,FITING)
        CALL FITST1(
     >              NUMKLE)
        NL2 = NUMKLE-1   ! FIT keyword is at position NUMKLE
        WRITE(6,201)
        CALL FITNU3(
     >              FITFNL)
C             write(*,*) ' main fitfnl, fiting : ',fitfnl, fiting 
        IF(FITFNL) THEN
          WRITE(6,200) 
          IF(NRES.GT.0) THEN
            WRITE(NRES,201)
            WRITE(NRES,200) 
 200        FORMAT(/,10X,
     >     ' MAIN PROGRAM :  FIT completed. ',
     >     ' Now doing final run using FIT variable values. ',A10)
          ENDIF
          ENDFIT = .FALSE.
          CALL FITNU8(FITFNL)
          CALL ZGOUBI(NL1,NL2,READAT,
     >                               NBLMN,ENDFIT)
          IF(NRES.GT.0) THEN 
            WRITE(NRES,201)
 201        FORMAT(/,132('*'))
            WRITE(NRES,334) NUMKLE,' Keyword FIT is skipped since '
     >      //'this is the (end of) final run following the fitting '
     >      //'procedure.','Now carrying on beyond FIT keyword.'
 334        FORMAT(/,2X,I5,2X,A,//,10X,A)
            CALL FLUSH2(NRES,.FALSE.)
          ENDIF
        ELSE
            WRITE(NRES,335) ' Final run following FIT[2] is skipped,'
     >      //'as requested.  Now carrying on beyond FIT keyword.'
 335        FORMAT(/,2X,A)
        ENDIF
        WRITE(6,201)
C Proceeds downstream of FIT[2] to the end of zgoubi.dat list
        READAT = .TRUE.
        FITING = .FALSE.
        CALL FITSTA(I6,FITING)
        NOEL = NUMKLE
        NL1 = NUMKLE + 1
        CALL ZGKLE(IQ(NL1-1),
     >                     KLE)   ! KLE = FIT !!
        CALL GO2KEY(NL1)
        NL2 = NBLMN
        CALL REBEL6(NL1, NBLMN)
        ENDFIT = .TRUE.
        CALL ZGOUBI(NL1,NL2,READAT,
     >                             NBLMN,ENDFIT)
        IF(.NOT. ENDFIT) THEN
          REWIND(NDAT)
          GOTO 11 
        ENDIF
      ENDIF

      GOTO 10
 
 996  WRITE(6,*) ' PGM ZGOUBI : error open file ', 
     >FILIN(DEBSTR(FILIN):FINSTR(FILIN))
      GOTO 10
 997  WRITE(6,*) ' PGM ZGOUBI : error open file ', 
     >FILOU(DEBSTR(FILOU):FINSTR(FILOU))
      GOTO 10
 995  WRITE(6,*) ' PGM ZGOUBI : error open file ', 
     >FILOG(DEBSTR(FILOG):FINSTR(FILOG))
      GOTO 10 

 99     CONTINUE
        WRITE(6   ,103) DMY,HMS
        IF(NRES.GT.0) THEN
          WRITE(NRES,103) DMY,HMS
          WRITE(NRES,101)
 101      FORMAT(/,128('*'))
          WRITE(NRES,FMT='(/,10X,
     >          '' Main program : stopped upon key  FIT'')')
        ENDIF

 10   CONTINUE
      
      IF(NRES.GT.0) THEN
        WRITE(NRES,fmt='(A)')  '   '
        WRITE(NRES,fmt='(A)')  '            Zgoubi run completed. '
        WRITE(NRES,103) DMY,HMS
      ENDIF

      WRITE(6   ,103) DMY,HMS
      CALL DATE2(DMY)
      CALL TIME2(HMS)
      IF(NRES.GT.0) WRITE(NRES,107) DMY,HMS
      WRITE(6   ,107) DMY,HMS
 107  FORMAT('  Job  ended  on    ',A,',  at  ',A,/)

      TEMP = TIMSEC
      CALL CPU_TIME(TIMSEC)
      IF(NRES.GT.0) WRITE(NRES,*) '  CPU time, total :  ',  TIMSEC-TEMP
      WRITE(   6,*) '  CPU time, total :  ',  TIMSEC-TEMP

      STOP
      END
