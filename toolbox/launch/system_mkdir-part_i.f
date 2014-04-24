      character cmmnd*200, txt400*400, txt9*9
      logical strcon
      integer debstr, finstr
      CHARACTER * 9   DMY,HMS

      CALL DATE2(DMY)
      CALL TIME2(HMS)

      cmmnd = 'ls -l | cat > echo'
      write(*,*) 'Pgm system_mkdir-part_i. Now executing '''
     >            ,cmmnd(debstr(cmmnd):finstr(cmmnd)),''''
      call system(cmmnd)
      
      open(unit=7,file = 'echo')
      
      ip = 1
 1    continue
        read(7,fmt='(a)',err=10,end=10) txt400
        txt400 = txt400(debstr(txt400):finstr(txt400))
        if(strcon(txt400,'part',4,
     >                              IS)) then 
!             write(*,*) ' is = ', is
          if(strcon(txt400(is+4:finstr(txt400)),'_',1,
     >                                                IS2)) then 
!             write(*,*) ' is, is2 = ', is, is2
!             write(*,*) 'txt400 = ',txt400(is+4:is+4+is2-2)
            read(txt400(is+4:is+4+is2-2),*,err=1,end=1) ipp
            if(ipp.gt.ip) ip = ipp
          endif                
        endif                
      goto 1

 10   continue
!             write(*,*) 'ip = ',ip 
        write(txt9,fmt='(i6)') ip+1

        cmmnd = 'mkdir -p '
     >  //'part'//txt9(debstr(txt9):finstr(txt9))//'_'
     >          //dmy//'.'//hms(debstr(hms):finstr(hms))
        write(*,*) 'Pgm system_mkdir-part_i. Now executing '''
     >            ,cmmnd(debstr(cmmnd):finstr(cmmnd)),''''
        call system(cmmnd)

        cmmnd = 'cp template.dat '
     >  //'part'//txt9(debstr(txt9):finstr(txt9))//'_'
     >          //dmy//'.'//hms(debstr(hms):finstr(hms))
     >  //'/zgoubi.dat'
        write(*,*) 'Pgm system_mkdir-part_i. Now executing '''
     >            ,cmmnd(debstr(cmmnd):finstr(cmmnd)),''''
        call system(cmmnd)

        cmmnd = '(cd '
     >  //'part'//txt9(debstr(txt9):finstr(txt9))//'_'
     >          //dmy//'.'//hms(debstr(hms):finstr(hms))
     >  //' ; ~/zgoubi/SVN/current/zgoubi/zgoubi '
     >  //' < zgoubi.dat > echo '
     >  //')'
        write(*,*) cmmnd
        call system(cmmnd)

      stop

      stop 
      end
      FUNCTION STRCON(STR,STRIN,NCHAR,
     >                                IS)
      implicit double precision (a-h,o-z)
      LOGICAL STRCON
      CHARACTER STR*(*), STRIN*(*)
C     ------------------------------------------------------------------------
C     .TRUE. if the string STR contains the string STRIN with NCHAR characters
C     at least once.
C     IS = position of first occurence of STRIN in STR
C     ------------------------------------------------------------------------

      INTEGER DEBSTR,FINSTR

      II = 0
      DO 1 I = DEBSTR(STR), FINSTR(STR)
        II = II+1
        IF( STR(I:I+NCHAR-1) .EQ. STRIN ) THEN
          IS = II
          STRCON = .TRUE.
          RETURN
        ENDIF
 1    CONTINUE
      STRCON = .FALSE.
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
      FUNCTION FINSTR(STRING)
      implicit double precision (a-h,o-z)
      INTEGER FINSTR
      CHARACTER * (*) STRING
C     --------------------------------------
C     RENVOIE DANS FINSTR LE RANG DU
C     DERNIER CHARACTER NON BLANC DE STRING,
C     OU BIEN 0 SI STRING EST VIDE ou BLANC.
C     --------------------------------------

      FINSTR=LEN(STRING)+1
1     CONTINUE
        FINSTR=FINSTR-1
        IF(FINSTR .EQ. 0) RETURN
        IF (STRING(FINSTR:FINSTR) .EQ. ' ') GOTO 1

      RETURN
      END
      include 'date2.f'
      include 'time2.f'










