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
      PROGRAM ZGOUBI_MAIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER * 9   DMY,HMS
      LOGICAL IDLUNI, READAT, FITING, ENDFIT
      CHARACTER(80) FNAME
      CHARACTER(80) FILEROOT
      integer ierri, filerootlength

      PARAMETER (I5=5, I6=6)

      PARAMETER (KSIZ=10)
      CHARACTER*(KSIZ) KLE
      integer debstr, finstr

      DATA ENDFIT / .FALSE. /

      CALL GET_COMMAND_ARGUMENT(1, FILEROOT, filerootlength, ierr)
      if(ierr .gt. 0) Then
        FILEROOT = "zgoubi"
        filerootlength=6
      endif
      
      if(filerootlength.gt.4) then
        if(FILEROOT(filerootlength-3:filerootlength) .eq. '.dat') then
          FILEROOT = FILEROOT(1:filerootlength-4)
          filerootlength = filerootlength -4 
        endif
      endif
      CALL INIDAT
      IF(IDLUNI(
     >          NDAT)) THEN
        FNAME = FILEROOT(1:filerootlength)//'.dat'
        OPEN(UNIT=NDAT,FILE=FNAME,STATUS='OLD',ERR=996)
c          write(*,*) ' zgoubi.dat is unit # ',ndat
      ELSE
        GOTO 996
      ENDIF

      IF(IDLUNI(
     >          NRES)) THEN
        FNAME = FILEROOT(1:filerootlength)//'.res'
        OPEN(UNIT=NRES,FILE=FNAME,ERR=997)
c          write(*,*) ' zgoubi.res is unit # ',nres
      ELSE
        GOTO 997
      ENDIF

      IF(IDLUNI(
     >          NLOG)) THEN
        FNAME = FILEROOT(1:filerootlength)//'.log'
        OPEN(UNIT=NLOG,FILE=FNAME,ERR=995)
      ELSE
        GOTO 995
      ENDIF

      CALL DATE2(DMY)
      CALL TIME2(HMS)
      CALL CPU_TIME(TIMSEC)

      WRITE(6   ,103) DMY,HMS
 103  FORMAT(/,'  Zgoubi, version 5.1.0.',/,
     >       '  Job  started  on  ',A,',  at  ',A)

 1    CONTINUE

      CALL INIDAT
      CALL RESET
      CALL CHECKS

      READAT = .TRUE.
      FITING = .FALSE.
      CALL FITSTA(I6,FITING)
      NL1 = 1
      NL2 = MXL
      ENDFIT = .FALSE.
      CALL ZGOUBI(NL1,NL2,READAT,NBEL,ENDFIT)

      CALL FITSTA(I5,
     >               FITING)
      IF(FITING) THEN
        READAT = .FALSE.
        CALL FITNU(NRES,*99)
        FITING = .FALSE.
        CALL FITSTA(I6,FITING)
        CALL FITST1(
     >              NUMKLE)
        NL2 = NUMKLE-1   ! FIT keyword is at position NUMKLE
        WRITE(6,201)
        WRITE(6,200) 
        IF(NRES.GT.0) THEN
          WRITE(NRES,201)
          WRITE(NRES,200) 
 200      FORMAT(/,10X,
     >   ' MAIN PROGRAM :  now final run using FIT values ',A10)
        ENDIF
        ENDFIT = .FALSE.
        CALL ZGOUBI(NL1,NL2,READAT,NBEL,ENDFIT)
c        write(*,*) ' zgoubi_main 2 fiting :',nl1,nl2,numkle
c        write(*,*) 
        WRITE(6,201)
        IF(NRES.GT.0) WRITE(NRES,201)
 201    FORMAT(/,128('*'))

C Proceeds until the end of zgoubi.dat list
        READAT = .TRUE.
        FITING = .FALSE.
        CALL FITSTA(I6,FITING)
        NOEL = NUMKLE
        NL1 = NUMKLE + 1
C        call go2key('FIT')
        CALL ZGKLE(IQ(NL1-1),
     >                     KLE)   ! KLE = FIT !!
        call go2key(NL1)
        NL2 = NBEL
        CALL REBEL6(NL1, NBEL)
        ENDFIT = .TRUE.
        CALL ZGOUBI(NL1,NL2,READAT,NBEL,ENDFIT)
c        IF(IDLUNI(
c     >          LUN)) THEN
c          OPEN(UNIT=LUN,FILE='zgoubi.fitSave',ERR=462)
c          CALL IMPAJU(LUN,F)          
c        ELSE
c          GOTO 462
c        ENDIF
c 462    CONTINUE

C For use of REBELOTE (REBELOTE encompasses FIT)
        IF(.NOT. ENDFIT) THEN
          REWIND(NDAT)
          GOTO 1 
        ENDIF
      ENDIF

      GOTO 10
 
 995  WRITE(6,*) ' PGM ZGOUBI : error open file ', FNAME
      GOTO 10
 996  WRITE(6,*) ' PGM ZGOUBI : error open file ', FNAME
      GOTO 10
 997  WRITE(6,*) ' PGM ZGOUBI : error open file ', FNAME
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

