      logical exs, idluni
      dimension borne(6), nbin(3)
      character txt*1

      data  kpa, kpb, kpstp / 1, 20000, 400/ 

      INQUIRE(FILE='lipsFitFromFai_iterate.In',exist=EXS)

      IF (IDLUNI(NLU)) THEN
        open(unit=nlu,file='lipsFitFromFai_iterate.In')
      ELSE
        stop 'Pgm lipsFitFromFai_iterate : no idle unit number ! '
      ENDIF

      if(.not. exs) then
       write(87,*)'WARNING : lipsFitFromFai_iterate.In does not exist'
       write(87,*)'Will be created from default values'

       write(nlu,*) kpa, kpb, kpstp, ' ! kpa, kpb, kpstp' 
      endif

      rewind(nlu)

        read(nlu,*,err=10,end=10) kpa, kpb, kpstp
        close(nlu)

C      call system('rm -f temp_ipmx')
      call system('mv -f lipsFitFromFai.Out  lipsFitFromFai.Out_old')
      call system('mv -f lipsFitFromFai_iterate.Out 
     >                             lipsFitFromFai_iterate.Out_old')

      do kp = kpa, kpb, kpstp
        IF (IDLUNI(NLU)) THEN
          open(unit=nlu,file='lipsFitFromFai.In')
        ELSE
          stop 'Pgm lipsFitFromFai_iterate :  No idle unit number - 2 !'
        ENDIF

        write(nlu,*) kp
        close(nlu)

        write(*,*) 
        write(*,*) ' ----------------------------------------------'
        write(*,*) 
     >  ' PGM lipsFitFromFai_iterate ; NOW NEW CALL TO lipsFitFromFai'
        write(*,*) ' Turn # : ',kp
        write(*,*) ' ----------------------------------------------'
        write(*,*) 
      call system('~/zgoubi/struct/tools/lipsFitFromFai/lipsFitFromFai')
      call system('cat lipsFitFromFai.Out >>lipsFitFromFai_iterate.Out')
        
c        IF (IDLUNI(itmp))  open(unit=itmp,file='temp_ipmx')
c        read(itmp,*) ierr,mxPss
c        close(itmp)
C        if(ierr .eq. -1) goto 11
      enddo
      close(nlu)

      stop

 10   stop 'Error during read in lipsFitFromFai_iterate.In'
C 11   write(*,*) 'End of .fai file reached, NPASS = ',mxPss
 11   write(*,*) 'End of .fai file reached '
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
