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
      use pariz_namelist_interface, only : MXX, MXY, IZ,
     >  initialize_input_parameters
      use xyzhc_interface, only : ensure_xyzhc_allocation
      use c_ss1_interface, only : ensure_xyz_allocation
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) LABEL
      INCLUDE "C.LABEL.H"     ! COMMON/LABEL/ LABEL(MXL,2)
      CHARACTER(10) DMY
      CHARACTER(9) HMS
      LOGICAL IDLUNI, READAT, FITING, FITBYD, FITRBL, FITNHB
      PARAMETER (I100 = 100)
      CHARACTER(I100) FLIN, FLOUT, FLOG

      CHARACTER(I100) FDAT
      SAVE FDAT

      PARAMETER (I5=5, I6=6)

      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLE
      LOGICAL FITFNL, FITLST
      INTEGER DEBSTR, FINSTR

      CHARACTER(LEN=I100), DIMENSION(:), ALLOCATABLE :: ARGS

      LOGICAL SAVXEC, SAVZPP
      LOGICAL OKWDAT, OKW

      CHARACTER(KSIZ) KEY, LBL1, LBL2
      PARAMETER (I0=0)

C Tentative pour faire fonctionner REBELOTE en compil dyn.
C Yet, still does not work. Try to remove one by one once REBELOTE fonctionne.
C (/home/meot/zgoubi/SVN/zgoubi-code/exemples/usersGuide/FIT-and-REBELOTE)
      SAVE FITING, READAT, NL1, NL2, FITBYD, NBLMN, NUMKLE,
     > FITFNL, FITRBL, OKWDAT, OKW, KLE, IPASS,NRBLT

C Dummy
      CHARACTER(1) TAB(1)
      SAVE IRET

      DATA FLIN, FLOUT, FLOG / 'zgoubi.dat', 'zgoubi.res', 'zgoubi.log'/
      DATA FDAT / 'zgoubi.dat' /
      DATA SAVXEC, SAVZPP / .FALSE., .FALSE.  /
      DATA TIMSEC / 0.D0 /
      DATA OKWDAT, OKW / .FALSE., .FALSE. /
      DATA TAB /  1 * ' '  /
      DATA FITFNL, FITLST / .FALSE., .FALSE. /
      DATA FITNHB / .FALSE. /                    ! .T. if FIt inhibited (by TWISS for instance)
      DATA IRANK  / 0 /
      DATA IRET / 0 /

      call initialize_input_parameters('pariz.nml')
      call ensure_xyzhc_allocation(MXX,MXY,IZ)
      call ensure_xyz_allocation(MXX,MXY,IZ)

C Manage possible arguments to zgoubi -----------------------
      NBARGS = COMMAND_ARGUMENT_COUNT()
      ALLOCATE(ARGS(NBARGS))
      DO IX = 1, NBARGS
         CALL GET_COMMAND_ARGUMENT(IX,ARGS(IX))
      END DO
      IF(NBARGS.GE.1) THEN
        WRITE(6,*) ' '
        WRITE(6,FMT='(A)') ' Argument(s) to zgoubi command :'
        DO IX = 1, NBARGS
          WRITE(6,FMT='(2(I0,A))') IX,' : '//ARGS(IX)
        END DO
      ENDIF
      IX = 1
      DO WHILE (IX .LE. NBARGS)
         IF(ARGS(IX) .EQ. '-rank'  ) THEN
           IX = IX + 1
           READ (ARGS(IX),'(I10)') IRANK
           CALL MCOBJB(IRANK)
         ELSEIF(  ARGS(IX) .EQ. '-fileIn'
     >   .OR.     ARGS(IX) .EQ. '-in'  ) THEN
           IX = IX + 1
           FLIN = ARGS(IX)
         ELSEIF(ARGS(IX) .EQ. '-fileOut'
     >   .OR.   ARGS(IX) .EQ. '-out'  ) THEN
           IX = IX + 1
           FLOUT = ARGS(IX)
         ELSEIF(ARGS(IX) .EQ. '-fileLog' ) THEN
           IX = IX + 1
           FLOG = ARGS(IX)
         ELSEIF(ARGS(IX) .EQ. '-saveExec') THEN
           SAVXEC = .TRUE.
         ELSEIF(ARGS(IX) .EQ. '-saveZpop') THEN
           SAVZPP = .TRUE.
         ELSE
           WRITE(*,*) 'argument # IX = ',IX,'/',NBARGS,' : ',ARGS(IX)
           WRITE(*,*) 'zgoubi_main : no such argument, # IX : ',ARGS(IX)
           CALL ENDJOB('Exiting.',-99)
         ENDIF
         IX = IX + 1
      ENDDO
C -----

      IF(FLIN .EQ. 'zgoubi_temp.dat') CALL ENDJOB(
     >'Pgm zgoubi_main. ''zgoubi_temp.dat'' input file name '//
     >'is reserved. Please choose a different one.')

      IF(IDLUNI(
     >          NLIN)) THEN
        FLIN = FLIN
        OPEN(UNIT=NLIN,FILE=FLIN,STATUS='OLD',ERR=996)
      ELSE
        GOTO 996
      ENDIF

      IF(IDLUNI(
     >          NRES)) THEN
      block 
        character(len=9) :: image_number
        character(len=:), allocatable :: output_file

        if (this_image()==1) then
          OPEN(UNIT=NRES,FILE=FLOUT,ERR=997)
        else
          write(image_number,'(i4)') this_image()
          output_file = trim(adjustl(FLOUT))
     >                       //"_image_"// trim(adjustl(image_number))
          OPEN(UNIT=NRES,FILE=output_file,ERR=997)
        end if
      end block
      ELSE
        GOTO 997
      ENDIF

      IF(IDLUNI(
     >          NLOG)) THEN
        OPEN(UNIT=NLOG,FILE=FLOG,ERR=995)
      ELSE
        GOTO 995
      ENDIF

      IF(SAVXEC .OR. SAVZPP) THEN
        WRITE(6,*) ' '
        WRITE(6,*) ' --'
        IF(SAVXEC) THEN
          WRITE(6,FMT='(A)')
     >    'Now copying zgoubi executable to local directory.'
          CALL SYSTEM('cp ~/zgoubi/SVN/current/zgoubi/zgoubi .')
        ENDIF
        IF(SAVZPP) THEN
          WRITE(6,FMT='(a)') 'Now copying zpop executable to local '
     >    //'directory (if you wish to run it, use an xterm terminal).'
          CALL SYSTEM('cp ~/zgoubi/SVN/current/zpop/zpop .')
        ENDIF
        WRITE(6,*) ' --'
        WRITE(6,*) ' '
      ENDIF

      CALL DATE2(DMY)
      CALL TIME2(HMS)
      CALL CPU_TIME(TIMSEC)

      WRITE(6   ,103) DMY,HMS
 103  FORMAT(/,'  Zgoubi, author''s dvlpmnt version.',/,
     >       '  Job  started  on  ',A,',  at  ',A)

      CALL PRDATA(NLIN,FLIN,FDAT,I100,NRES,
     >                                LABEL,NBLM,NDAT)
      NBLMN = NBLM

      READAT = .TRUE.

 11   CONTINUE
      FITING = .FALSE.
      CALL INIDAT
      CALL RESET
      CALL CHECKS
      CALL FITSTA(I6,FITING)
      NL1 = 1
      NL2 = MXL
      FITBYD = .FALSE.
      CALL FITST4(FITBYD)
      NBLMI = NBLMN
C FM - 17.10.24. Allows continuing beyond FIT
 12   CONTINUE
      CALL ZGOUBI(NL1,NL2,READAT,NBLMI)
      CALL ZGIRET(
     >            IRET)
      IF(IRET.EQ.1) GOTO 10
      CALL FITSTA(I5,
     >               FITING)
      IF(.NOT. FITNHB) THEN
       IF(FITING) THEN
        IF(NRES.GT.0) WRITE(NRES,FMT='(5X,
     >  ''FIT procedure launched.'',/)')
C     >  ''FIT procedure launched. Method is '',I1,/)') MTHOD
        READAT = .FALSE.
        CALL FITNU(NRES)
        FITING = .FALSE.
        CALL FITSTA(I6,FITING)
        CALL FITST1(
     >              NUMKLE)
        NL2 = NUMKLE-1   ! FIT keyword is at position NUMKLE
        WRITE(6,201)
        CALL FITST5(
     >              FITFNL)
        IF(FITFNL) THEN
          FITLST = .TRUE.
          CALL FITSTX(FITLST)
          WRITE(6,200)
          IF(NRES.GT.0) THEN
            WRITE(NRES,201)
            WRITE(NRES,200)
 200        FORMAT(/,10X,
     >     ' MAIN PROGRAM :  FIT completed. ',
     >     ' Now doing a last run using variable values from FIT. ',A10)
          ENDIF
          FITBYD = .FALSE.
          CALL FITST4(FITBYD)
          NBLMI = NBLMN
          CALL ZGOUBI(NL1,NL2,READAT,NBLMI)
          FITLST = .FALSE.
          CALL FITSTX(FITLST)
          IF(NRES.GT.0) THEN
            WRITE(NRES,201)
 201        FORMAT(/,132('*'))
            WRITE(NRES,334) NUMKLE,' Keyword FIT[2] is skipped since '
     >      //'this is the (end of) last run following the fitting '
     >      //'procedure.','Now carrying on beyond FIT keyword.'
 334        FORMAT(/,2X,I5,2X,A,//,10X,A,/)
            CALL FLUSH2(NRES,.FALSE.)
          ENDIF

        ELSE
          IF(NRES.GT.0)
     >    WRITE(NRES,335) ' Last run following FIT[2] is skipped,'
     >    //' as requested.  Now carrying on beyond FIT keyword.'
 335      FORMAT(/,2X,A)
        ENDIF

        OKWDAT = .TRUE.
        OKW = OKWDAT

        WRITE(6,201)
C Proceeds downstream of FIT[2]  toward end of zgoubi.dat list (possibly meeting REBELOTE)
        READAT = .TRUE.
        FITING = .FALSE.
        CALL FITSTA(I6,FITING)
        NOEL = NUMKLE
        NL1 = NUMKLE + 1
        CALL ZGKLE(IQ(NL1-1),
     >                     KLE)   ! KLE = FIT !!
        CALL GO2KEY(NL1,I0,I0,TAB,
     >                            KEY,LBL1,LBL2)
        NL2 = NBLMN
        NBLMI = NBLMN
C        CALL REBEL6(NL1, NBLMI)
        CALL REBLT6(NL1, NBLMI)
        FITBYD = .TRUE.
        CALL FITST4(FITBYD)
        NBLMI = NBLMN
        CALL ZGOUBI(NL1,NL2,READAT,NBLMI)
        CALL ZGIRET(
     >            IRET)
        IF(IRET.EQ.1) GOTO 10
        CALL FITST7(
     >              FITRBL)   ! Switched to T by REBELOTE, SVDOC... if FIT embedded
        IF(FITRBL) THEN
          OKW = OKWDAT
          IF(OKWDAT) THEN
            CALL FITWDA(
     >                  IER)
            OKWDAT = .FALSE.
          ENDIF
          CALL ZGIPAS(
     >                IPASS,NRBLT)
          IF(IPASS .LE. NRBLT+1) THEN
            READAT = .FALSE.
            FITRBL = .FALSE.
            CALL FITST8(FITRBL)
            REWIND(NDAT)
            GOTO 11
          ELSE
            GOTO 10
          ENDIF
        ENDIF

       ENDIF
      ENDIF

C FM - 17.10.24. FIT has possibly been inhibited (e.g., by TWIIS)
      CALL FITSTB(
     >            FITNHB)
C FM - 17.10.24. Allows carrying on beyond FIT
      GOTO 12
C      GOTO 10

 996  WRITE(6,*) ' PGM ZGOUBI : error open file ',
     >FLIN(DEBSTR(FLIN):FINSTR(FLIN))
      GOTO 10
 997  WRITE(6,*) ' PGM ZGOUBI : error open file ',
     >FLOUT(DEBSTR(FLOUT):FINSTR(FLOUT))
      GOTO 10
 995  WRITE(6,*) ' PGM ZGOUBI : error open file ',
     >FLOG(DEBSTR(FLOG):FINSTR(FLOG))
      GOTO 10

  10   CONTINUE

      CLOSE(NDAT)

      IF(NRES.GT.0) THEN
        WRITE(NRES,FMT='(A)')  '   '
        WRITE(NRES,FMT='(A)')  '            ZGOUBI RUN COMPLETED. '
        WRITE(NRES,103) DMY,HMS
      ENDIF

      WRITE(6   ,103) DMY,HMS
      CALL DATE2(DMY)
      CALL TIME2(HMS)
      IF(NRES.GT.0) WRITE(NRES,107) DMY,HMS
      WRITE(6   ,107) DMY,HMS
 107  FORMAT('  JOB  ENDED  ON    ',A,',  AT  ',A,/)

      TEMP = TIMSEC
      CALL CPU_TIME(TIMSEC)
      IF(NRES.GT.0) WRITE(NRES,*) '  CPU time, total :  ',  TIMSEC-TEMP
      WRITE(   6,*) '  CPU time, total :  ',  TIMSEC-TEMP

      IF(OKW) CALL FITWDA(
     >                    IER)

      CLOSE(ABS(NRES))
      CLOSE(NLOG)

C           CALL ENDJOB
C     >     ('Pgm zgoubi_main : Job ended.',-9999)

      END
