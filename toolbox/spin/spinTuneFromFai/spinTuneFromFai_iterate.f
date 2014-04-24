      logical exs, idluni

      data  kpa, kpb, kpstp, ksmpl / 1, 60000, 1000, 200/ 
      data born1, born2 / 0.4, .6 /
      data nbin / 2500 /


      INQUIRE(FILE='spinTuneFromFai_iterate.In',exist=EXS)

      IF (IDLUNI(NLU)) THEN
          open(unit=nlu,file='spinTuneFromFai_iterate.In')
      ELSE
          stop 'Pgm spinTuneFromFai_iterate : no idle unit number !'
      ENDIF

      if(.not. exs) then
       write(*,*)
     >  'WARNING: File spinTuneFromFai_iterate.In does not exist'
       write(*,*)
     >  'spinTunesFromFai_iterate creates one from default values'

       write(nlu,*) kpa, kpb, kpstp, ksmpl, ' ! kpa, kpb, kpstp, ksmpl' 
       write(nlu,*) born1, born2, ' ! borne 1, borne 2' 
       write(nlu,*) nbin, ' ! nbin '
      endif

      rewind(nlu)

        read(nlu,*,err=10,end=10) kpa, kpb, kpstp, ksmpl
        read(nlu,*,err=10,end=10) born1, born2
        read(nlu,*,err=10,end=10) nbin
        close(nlu)

      call system('cat spinTuneFromFai.Out >>spinTuneFromFai.Out_old')
      call system('cat spinTuneFromFai_spectrum.Out 
     >                            >>spinTuneFromFai_spectrum.Out_old')
      call system('cat spinTuneFromFai_iterate.Out >> 
     >                           spinTuneFromFai_iterate.Out_old')
      call system('cat spinTuneFromFai_spectrum_iterate.Out >> 
     >                      spinTuneFromFai_spectrum_iterate.Out_old')
      call system('rm -f spinTuneFromFai.Out ')
      call system('rm -f spinTuneFromFai_spectrum.Out ')
      call system('rm -f spinTuneFromFai_iterate.Out')
      call system('rm -f spinTuneFromFai_spectrum_iterate.Out')

      do kp = kpa, kpb, kpstp
        IF (IDLUNI(NLU)) THEN
          open(unit=nlu,file='spinTuneFromFai.In')
        ELSE
          stop 'Pgm spinTuneFromFai_iterate :   No idle unit number ! '
        ENDIF

        write(nlu,fmt='(2(i6,1x),t60,t60,a)')  kp, ksmpl
     >  ,' ! kpa, ksmpl : Fourier transf. ksmpl turns starting from kpa'
        write(nlu,fmt='(a,1x,t60,a)')  ' 1' 
     >  ,' ! itraj : number of the particle to be Fourier''ed'
        write(nlu,fmt='(2(f8.4,1x),t60,a)')  born1, born2
     >  ,' ! Q_1, Q_2 : x/y/l  spectrum range'
        write(nlu,fmt='(i6,1x,t60,a)')  nbin
     >  ,' ! nbin : x/y/l # of bins in spectrum range'
        close(nlu)

        write(*,*) 
        write(*,*) ' ----------------------------------------------'
        write(*,*) 
     >  ' PGM spinTuneFromFai_iterate ; NOW NEW CALL spinTuneFromFai'
        write(*,*) ' Turn #, range : ',kp,'-',kp+ksmpl-1
        write(*,*) ' ----------------------------------------------'
        write(*,*) 
        call system
     >  ('~/zgoubi/struct/tools/spin/spinTuneFromFai/spinTuneFromFai')
        call system('cat spinTuneFromFai.Out 
     >                       >> spinTuneFromFai_iterate.Out')
        call system('cat spinTuneFromFai_spectrum.Out 
     >                       >> spinTuneFromFai_spectrum_iterate.Out')
      enddo
      close(nlu)

      stop

 10   stop 'Error during read in spinTuneFromFai_iterate.In'
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
