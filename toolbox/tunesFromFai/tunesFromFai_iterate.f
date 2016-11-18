      logical exs, idluni, oksav
      dimension borne(6), nbin(3)
      character txt*1,  txtQ
      logical okQ

      character(200) TXT200
      CHARACTER(80) STRA(5)
      CHARACTER(40) status

      logical strcon

C RHIC
C      data  kpa, kpb, kpstp, ksmpl / 1, 80000, 2000, 200/ 
C      data borne / 0.5, 1., 0.5, 1., 0.5, 1. /
C AGD
      data  kpa, kpb, kpstp, ksmpl / 1, 10000, 400, 400/ 
      data borne / 0.5, 1., 0.5, 1., 0., 1. /
      data nbin / 3*2500 /
      data oksav / .false. /
      data txt / 'n' /
      data txtQ / 'y' /
      data okQ / .true. /
      data nt / -1 /
      data stra / 5*'  ' /

      INQUIRE(FILE='tunesFromFai_iterate.In',exist=EXS)

      IF (IDLUNI(NLU)) THEN
        open(unit=nlu,file='tunesFromFai_iterate.In')
      ELSE
        stop 'Pgm tunesFromFai_iterate : no idle unit number ! '
      ENDIF

      if(.not. exs) then
       write(87,*)'WARNING: File tunesFromFai_iterate.In does not exist'
       write(87,*)'tunesFromFai_iterate creates one from default values'

       write(nlu,*) kpa, kpb, kpstp, ksmpl, nt,
     > ' ! kpa, kpb, kpstp, ksmpl, nt ' 
       write(nlu,*) (borne(i),i=1,6), ' ! (borne 1, borne 2) * 3' 
       write(nlu,*) (nbin(i),i=1,3), ' ! nbin '
       write(nlu,*) txt,' !  save spectra (y/n)'
       write(nlu,*) txtQ,' ! yes/no compute tunes, too (not just lips)'
      endif

      rewind(nlu)

C        read(nlu,*,err=10,end=10) kpa, kpb, kpstp, ksmpl
        read(nlu,fmt='(a)',err=10,end=10) txt200
        if(strcon(txt200,'!',
     >                       is)) txt200 = txt200(1:is-1)
        call strget(txt200,5,
     >                       nbstr,stra)
        read(stra(1),*) kpa
        read(stra(2),*) kpb
        read(stra(3),*) kpstp
        read(stra(4),*) ksmpl
        if(nbstr.eq.5) read(stra(5),*) nt
        read(nlu,*,err=10,end=10) (borne(i),i=1,6)
        read(nlu,*,err=10,end=10) (nbin(i),i=1,3)
C        read(nlu,fmt='(L1)',err=10,end=10) oksav
        read(nlu,*,err=10,end=10) txt
        read(nlu,*,err=33,end=33) txtQ
 33     if(txtQ .ne. 'n') txtQ = 'y'
        oksav = txt.eq. 'y'
        okQ = txtQ .eq. 'y'
        close(nlu)

      call system('rm -f temp_ipmx')
      call system('mv -f tunesFromFai.out  tunesFromFai.out_old')
      call system('mv -f tunesFromFai_spctra.Out 
     >                              tunesFromFai_spctra.Out_old')
      call system('mv -f tunesFromFai_iterate.out 
     >                             tunesFromFai_iterate.out_old')
      call system('mv -f tunesFromFai_iterate_spctra.Out  
     >                      tunesFromFai_iterate_spctra.Out_old')

      do kp = kpa, kpb, kpstp
        IF (IDLUNI(NLU)) THEN
          open(unit=nlu,file='tunesFromFai.In')
        ELSE
          stop 'Pgm tunesFromFai_iterate :   No idle unit number - 2 ! '
        ENDIF

        write(nlu,*) kp, ksmpl, nt
        write(nlu,*) (borne(i),i=1,6), ' ! (borne 1, borne 2) * 3' 
        write(nlu,*) (nbin(i),i=1,3), ' ! nbin '
C        write(nlu,*)  oksav,'    !  save spectra '
        write(nlu,*)  txt
        write(nlu,*) txtQ
        oksav = txt.eq. 'y'
        okQ = txtQ.eq. 'y'
        close(nlu)

        write(*,*) 
        write(*,*) ' ----------------------------------------------'
        write(*,*) 
     >  ' PGM tunesFromFai_iterate ; NOW NEW CALL TO tunesFromFai'
        write(*,*) ' Turn #, range : ',kp,'-',kp+ksmpl-1
        write(*,*) ' ----------------------------------------------'
        write(*,*) 
        call system
     >  ('~/zgoubi/SVN/current/toolbox/tunesFromFai/tunesFromFai')

        if(idluni(itmp)) open(unit=itmp,file='tunesFromFai.status')
        read(itmp,fmt='(a)') status
               write(*,*) ' iterate status : ',status
        if(strcon(status,'FAILED'
     >                           ,IS)) stop 
     >         'Pgm tunesFromFia_iterate. Cannot find .fai file.'
        close(itmp)        

        write(*,*) ' '
        write(*,*) 'Pgm lipsFitFromFai_iterate. Now doing :'
        write(*,*) 'cat tunesFromFai.out >> tunesFromFai_iterate.out'
        write(*,*) ' '
C           read(*,*)

        call system('cat tunesFromFai.out >> tunesFromFai_iterate.out')
        call system('cat tunesFromFai_spctra.Out >>
     >                            tunesFromFai_iterate_spctra.Out')
        
        IF (IDLUNI(itmp))  open(unit=itmp,file='temp_ipmx')
        read(itmp,*) ierr,mxPss
        close(itmp)
C        if(ierr .eq. -1) goto 11
C        if(mxPss .le. 0) goto 12

      enddo
      close(nlu)

      stop

 10   stop 'Pgm tunesFromFai_iterate. '
     >//'Error during read in tunesFromFai_iterate.In'
 11   write(*,*) 'Pgm tunesFromFai_iterate. '
      write(*,*) 'Error report from pgm tunesFromFai. '
      write(*,*) 'mxPss = ',mxPss
      stop
 12   write(*,*) 'Pgm tunesFromFai_iterate. '
      write(*,*) 'End of .fai file reached, found mxPss = ',mxPss
      stop
      end
      FUNCTION IDLUNI(LN)
      LOGICAL IDLUNI

      LOGICAL OPN

      I = 20
 1    CONTINUE
        INQUIRE(UNIT=I,ERR=99,IOSTAT=IOS,OPENED=OPN)
        I = I+1
        IF(I .EQ. 100) GOTO 99
        IF(OPN) GOTO 1
        IF(IOS .GT. 0) GOTO 1
      
      LN = I-1
      IDLUNI = .TRUE.
      RETURN

 99   CONTINUE
      LN = 0
      IDLUNI = .FALSE.
      RETURN
      END
      FUNCTION STRCON(STR,STR2,
     >                         IS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL STRCON
      CHARACTER STR*(*), STR2*(*)
C     ---------------------------------------------------------------
C     .TRUE. if the string STR contains the string STR2 at least once
C     IS = position of first occurence of STR2 in STR 
C     (i.e.,STR(IS:IS+LEN(STR2)-1)=STR2)
C     ---------------------------------------------------------------
      INTEGER DEBSTR,FINSTR
      LNG2 = LEN(STR2(DEBSTR(STR2):FINSTR(STR2)))
      IF(LEN(STR).LT.LNG2) GOTO 1
      DO I = DEBSTR(STR), FINSTR(STR)-LNG2+1
        IF( STR(I:I+LNG2-1) .EQ. STR2 ) THEN
          IS = I 
          STRCON = .TRUE.
          RETURN
        ENDIF
      ENDDO
 1    CONTINUE
      STRCON = .FALSE.
      RETURN
      END
      SUBROUTINE STRGET(STR,MSS,
     >                          NST,STRA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) STR, STRA(MSS)
C     ------------------------------------------------------
C     Extract substrings #1 up to #MSS, out of string STR. 
C     Strings are assumed spaced by (at least) one blank. 
C     They are saved in  array STRA, and their actual number 
C     (possibly < mss) is NST.
C     ------------------------------------------------------
      INTEGER FINSTR

      CHARACTER STR0*(300)

      if(len(str0).lt.len(str)) 
     >  stop ' SBR STRGET : Increase length of string str0'

      STR0 = STR
      IE = FINSTR(STR)
      NST = 0
      I2 = 1

 1    CONTINUE

        IF(STR(I2:I2) .EQ. ' '  .OR. 
     >     STR(I2:I2) .EQ. ',') THEN
          I2 = I2 + 1
          IF(I2 .LE. IE) GOTO 1
        ELSE
          I1 = I2
 2        CONTINUE
          I2 = I2 + 1
          IF(I2 .LE. IE) THEN
            IF(STR(I2:I2) .EQ. ' '  .OR. 
     >         STR(I2:I2) .EQ. ',') THEN
              IF(NST .LT. MSS) THEN
                NST = NST + 1
                STRA(NST) = STR(I1:I2-1)
                I2 = I2 + 1
                GOTO 1
              ENDIF
            ELSE
              GOTO 2
            ENDIF
          ELSE
            IF(STR(I2-1:I2-1) .NE. ' ' .AND.
     >         STR(I2-1:I2-1) .NE. ',') THEN
              IF(NST .LT. MSS) THEN
                NST = NST + 1
                STRA(NST) = STR(I1:I2-1)
              ENDIF
            ENDIF
          ENDIF
        ENDIF

      STR = STR0

      RETURN
      END
      FUNCTION DEBSTR(STRING)
      implicit double precision (a-h,o-z)
      INTEGER DEBSTR
      CHARACTER * (*) STRING

C     --------------------------------------
C     RENVOIE DANS DEBSTR LE RANG DU
C     1-ER CHARACTER NON BLANC DE STRING,
C     OU BIEN 0 SI STRING EST VIDE ou BLANC.
C     --------------------------------------

      DEBSTR=0
      LENGTH=LEN(STRING)
C      LENGTH=LEN(STRING)+1
1     CONTINUE
        DEBSTR=DEBSTR+1
C        IF(DEBSTR .EQ. LENGTH) RETURN
C        IF (STRING(DEBSTR:DEBSTR) .EQ. ' ') GOTO 1
        IF (STRING(DEBSTR:DEBSTR) .EQ. ' ') THEN
          IF(DEBSTR .EQ. LENGTH) THEN
            DEBSTR = 0
            RETURN
          ELSE
            GOTO 1
          ENDIF
        ENDIF

      RETURN
      END
      FUNCTION FINSTR(STR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER FINSTR
      CHARACTER * (*) STR
C     -----------------------------------
C     Renvoie dans FINSTR le rang du
C     dernier caractere non-blanc de STR.
C     Renvoie 0 si STR est vide ou blanc.
C     -----------------------------------

      FINSTR=LEN(STR)+1
1     CONTINUE
         FINSTR=FINSTR-1
         IF(FINSTR.EQ. 0) RETURN
         IF (STR(FINSTR:FINSTR).EQ. ' ') GOTO 1
      RETURN
      END
