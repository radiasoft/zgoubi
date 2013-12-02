      logical exs, idluni, oksav
      dimension borne(6), nbin(3)
      character txt*1,  txtQ
      logical okQ

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

      INQUIRE(FILE='tunesFromFai_iterate.In',exist=EXS)

      IF (IDLUNI(NLU)) THEN
        open(unit=nlu,file='tunesFromFai_iterate.In')
      ELSE
        stop 'Pgm tunesFromFai_iterate : no idle unit number ! '
      ENDIF

      if(.not. exs) then
       write(87,*)'WARNING: File tunesFromFai_iterate.In does not exist'
       write(87,*)'tunesFromFai_iterate creates one from default values'

       write(nlu,*) kpa, kpb, kpstp, ksmpl, ' ! kpa, kpb, kpstp, ksmpl' 
       write(nlu,*) (borne(i),i=1,6), ' ! (borne 1, borne 2) * 3' 
       write(nlu,*) (nbin(i),i=1,3), ' ! nbin '
       write(nlu,*) txt,' !  save spectra (y/n)'
       write(nlu,*) txtQ,' ! yes/no compute tunes, too (not just lips)'
      endif

      rewind(nlu)

        read(nlu,*,err=10,end=10) kpa, kpb, kpstp, ksmpl
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

        write(nlu,*) kp, ksmpl
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
        call system('cat tunesFromFai.out >> tunesFromFai_iterate.out')
        call system('cat tunesFromFai_spctra.Out >>
     >                            tunesFromFai_iterate_spctra.Out')
        
        IF (IDLUNI(itmp))  open(unit=itmp,file='temp_ipmx')
        read(itmp,*) ierr,mxPss
        close(itmp)
C        if(ierr .eq. -1) goto 11
      enddo
      close(nlu)

      stop

 10   stop 'Error during read in tunesFromFai_iterate.In'
 11   write(*,*) 'End of .fai file reached, NPASS = ',mxPss
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
