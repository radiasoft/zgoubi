      logical exs, idluni
      dimension borne(6), nbin(3)
      character(1) txt
      character(200) txt200
      character(80) filfai
      character(1) rmvDsp1, rmvDspi, rmvDsp
      logical strcon
      integer debstr, finstr

      data filfai / 'zgoubi.fai' /
      data  kpa, kpb, kpstp / 1, 9, 1/ 
      data  kla, klb, klstp / -1, -1, -1 / ! Will not consider lmnt #. Otherwise kpa<kpb. 
      data rmvDsp / 'N' /

      INQUIRE(FILE='lipsFitFromFai_iterate.In',exist=EXS)

      IF (IDLUNI(lIn)) THEN
        open(unit=lIn,file='lipsFitFromFai_iterate.In')
      ELSE
        stop 'Pgm lipsFitFromFai_iterate : no idle unit number ! '
      ENDIF

      if(.not. exs) then
       write(87,*)'WARNING : lipsFitFromFai_iterate.In does not exist'
       write(87,*)'Will be created from default values'

       write(lIn,*) filfai,          ' ! .fai file name' 
       write(lIn,*) kpa, kpb, kpstp, ' ! kpa, kpb, kpstp (pass #)' 
       write(lIn,*) kla, klb, klstp, ' ! kla, klb, klstp (lmnt #)' 
       write(lIn,fmt='(t11,a,1x,t30,a)')  rmvDsp
     >  ,' ! rmvDsp : remove dispersion Y/N (default is No), '
     >  //' Y requires zgoubi.TWISS.out from prior TWISS run.'
      endif

      rewind(lIn)

        read(lIn,*,err=10,end=10) filfai
        read(lIn,*,err=10,end=10) kpa, kpb, kpstp
        kla1 = kla ; klb1 = klb ; klstp1 = klstp
        read(lIn,*,err=20,end=20) klai, klbi, klstpi
        kla1 = klai ; klb1 = klbi ; klstp1 = klstpi
 20     continue
        kla = klai ; klb = klbi ; klstp = klstpi
        rmvDsp1 = rmvDsp
        read(lIn,*,err=21,end=21) txt200
        if(strcon(txt200,'!',
     >                       is)) txt200 = txt200(debstr(txt200):is-1)
        read(txt200(1:1),*,err=21,end=21) rmvDspi
        rmvDsp1 = rmvDspi
 21     continue
        rmvDsp = rmvDsp1
        if(rmvDsp .ne. 'Y') rmvDsp = 'N'

       rewind(lIn)
       write(lIn,fmt='(2a)') filfai,          ' ! .fai file name' 
       write(lIn,*) kpa, kpb, kpstp, ' ! kpa, kpb, kpstp (pass #)' 
       write(lIn,*) kla, klb, klstp, ' ! kla, klb, klstp (lmnt #)' 
       write(lIn,fmt='(t11,a,1x,t30,a)')  rmvDsp
     >  ,' ! rmvDsp : remove dispersion Y/N (default is No), '
     >  //' Y requires zgoubi.TWISS.out from prior TWISS run.'

        close(lIn)

C      call system('rm -f temp_ipmx')
      call system('mv -f lipsFitFromFai.Out  lipsFitFromFai.Out_old')
      call system('mv -f lipsFitFromFai_iterate.Out 
     >                             lipsFitFromFai_iterate.Out_old')

      do kp = kpa, kpb, kpstp
       do kl = kla, klb, klstp
        IF (IDLUNI(lIn)) THEN
          open(unit=lIn,file='lipsFitFromFai.In')
        ELSE
          stop 'Pgm lipsFitFromFai_iterate :  No idle unit number - 2 !'
        ENDIF

        write(lIn,*) filfai
        write(lIn,*) kp
        write(lIn,*) kl 
        write(lIn,*) rmvDsp
        close(lIn)

        write(*,*) 
        write(*,*) ' ----------------------------------------------'
        write(*,*) 
     >  ' PGM lipsFitFromFai_iterate ; NOW NEW CALL TO lipsFitFromFai'
        write(*,*) ' Turn # : ',kp
        write(*,*) ' Lmnt # : ',kl
        write(*,*) ' Dispersion removal : ',rmvDsp
        write(*,*) ' ----------------------------------------------'
        write(*,*) 
        call system
     >  ('~/zgoubi/SVN/current/toolbox/lipsFitFromFai/lipsFitFromFai')
        call system('cat lipsFitFromFai.Out '
     >                  //' >> lipsFitFromFai_iterate.Out')
        
c        IF (IDLUNI(itmp))  open(unit=itmp,file='temp_ipmx')
c        read(itmp,*) ierr,mxPss
c        close(itmp)
C        if(ierr .eq. -1) goto 11
       enddo
      enddo
      close(lIn)

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
