C Starting from initial coordinates near closed orbits, as read from input file 
C zgoubi_StabLim-In.dat (e.g., a copy of zgoubi_geneMap-Out.dat or zgoubi_searchCO-Out.dat), 
C will search stability limits. 
C zgoubi.dat is presumed to contain stable trajectories, which will be pushed to 
C extreme amplitude
      parameter (lunR=11,lunW=12,lunIn=15)
      character txt132*132, txt20*6, let*1, namFil*30
      parameter (nTrajmx=999)
      dimension x(nTrajmx),xp(nTrajmx),z(nTrajmx),
     >                   zp(nTrajmx),s(nTrajmx),d(nTrajmx),let(nTrajmx)
      dimension xi(nTrajmx),xpi(nTrajmx),zi(nTrajmx),
     >                   zpi(nTrajmx)
      dimension storb(6,nTrajmx)
      logical ok, result, strcon, first
      character*2 HV, HVIn
      character*2 txt0
      integer debstr, finstr
      logical  empty

      data lunData, icho / 7, 10 /
C H=pure H, Hz=H+small z, V=vertical 
      data txt0 / '0 ' /
      data HV / 'Hz' /
      data prec / 1.d-2 / !cm
      data duStrt / 0.5 /  ! cm
      data nTurn / 499 / 
      data kTrk / 0 / 
      data njj / 5 /
      data kprx / 1 / 
      data ok, first/ .false., .true. /
      data kobj2 / -999 /

      call system('mv -f searchStabLim-tunes.out 
     >                    searchStabLim-tunes.out_old')

      open(unit=lunData,file='searchStabLim.data')
      open(unit=icho,file='searchStabLim.echo')

 111  continue

C---------- Input data ------------------
      read(lunData,*,err=66,end=998) 
     >       HVIn, precIn, duStI, nTurnI, kTrkIn, njjIn, kprxI
c      close(lunData)
      HV = HVIn
      if(HV.ne.'H' .and. HV.ne.'Hz' .and. HV.ne.'V') goto 997
      prec = precIn  ! cm
      kTrk = kTrkIn  ! will track for tunes
      njj = njjIn    ! number of ellipses to be tracked after stab limit is found 
      duStrt = duStI ! = dx or dz starting increment
      nTurn = nTurnI
      kprx = kprxI   ! 1 to have paraxial traj tracked in addition to njj
 66   continue
      write(6,*) '  H/V, precision (cm), #turns  : ',HV, prec, nTurn
      write(icho,*) '  H/V, precision (cm), #turns  : ',HV, prec, nTurn
C               stop
      call flush2(icho,.false.)
C----------------------------------------
      result = .false.

      call system('mv -f zgoubi.dat zgoubi.dat-copy')

C zgoubi_StabLim-In.dat is supposed to contain stable trajectories 
C - these could be output from prior searchCO procedure - 
C these will be starting point for search of stability limit. 
C
      open(unit=lunR,file='zgoubi_StabLim-In.dat')


      jo = 1
      jok = 0
      xst = 0.
 2    continue
 
        open(unit=lunW,file='zgoubi.dat')
        rewind(lunW)
        rewind(lunR)

C Read till "KOBJ"
c        do i=1,3
c          read(lunR,fmt='(a)') txt132
c          write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c        enddo
        txt132 = ' '
        dowhile
     >  (txt132(2:6) .ne. 'OBJET' .and. TXT132(2:8) .ne. 'MCOBJET')
          read(lunR,fmt='(a)') txt132
          txt132 = txt132(debstr(txt132):finstr(txt132))
          if(empty(txt132)) then
            write(lunW,fmt='(a)') '!'
          else
             write(lunW,fmt='(a)')
     >        txt132(debstr(txt132):finstr(txt132))
          endif
c          write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
        enddo
        read(lunR,fmt='(a)') txt132
          if(empty(txt132)) then
            write(lunW,fmt='(a)') '!'
          else
             write(lunW,fmt='(a)')
     >        txt132(debstr(txt132):finstr(txt132))
          endif
c        write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
        read(lunR,fmt='(a)') txt132
        write(lunW,*) ' 2 '
C Set "IMAX IMAXT" for kobj=2 option 
        read(txt132,*) xobj
        kobj = int(xobj)
        if    (kobj .eq. 2) then
          read(lunR,*) nTraj, imaxt
        elseif(kobj .eq. 5) then
          txt132 = txt132(debstr(txt132):finstr(txt132))
          read(txt132(3:99),*) kobj2
          ntraj = kobj2
          read(lunR,fmt='(a)') txt132        ! skip sample
        endif

c               write(*,*) ' Pgm searchStabLim. ntraj,kobj,kobj2 : ',
c     >           ntraj,kobj,kobj2
c                    read(*,*)

        if(nTraj.gt.nTrajmx) 
     >  stop 'Pgm searchStabLim. Too many trajectories...'

        txt132 = ' 1  1'
        write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
C Read all initial traj from zgoubi_StabLim-In.dat, supposed to be stable
        do i=1,nTraj
          read(lunR,*) x(i),xp(i),z(i),zp(i),s(i),d(i),let(i)
          xi(i) =x(i); xpi(i) =xp(i); zi(i) =z(i); zpi(i) =zp(i)
        enddo

C Retains only one initial traj at a time
C and sets z to nul
        if    (HV .eq. 'H') then
          z(jo) = 0.
        elseif(HV .eq. 'Hz') then
          z(jo) = 1.e-4
        elseif(HV .eq. 'V') then
          z(jo) = 0.
        elseif(HV .eq. 'HV') then
          z(jo) = xst
        else
          write(6,*)
          write(6,*) 'Last value of HV read in searchStabLim.data is '''
     >                  , HV,''''
          write(6,*) 'End of job... '
          goto 99
        endif
        zp(jo) = 0.
        write(lunW,fmt='(1p,4e16.8,e9.1,e16.8,4a,i4)') 
     >    x(jo),xp(jo),z(jo),zp(jo),s(jo),d(jo),' ','''',let(jo),'''',jo
        write(6,*) 
        write(6,*) ' ---------    New trajectory :'
        write(6,fmt='(1p,4e16.8,e9.1,e16.8,4a,i4)') 
     >    x(jo),xp(jo),z(jo),zp(jo),s(jo),d(jo),' ','''',let(jo),'''',jo
        write(6,*)

C Complete OBJET with the line with '1'
        if(kobj .eq. 2) then
 37       read(lunR,fmt='(a)',end=62) txt132
          txt132 = txt132(debstr(txt132):finstr(txt132))
          if(txt132(1:2) .eq. '1 ') then
            goto 37
          else
            backspace(lunR)
          endif
        endif
        write(lunW,*) ' 1 '
        txt132 = '''FAISTORE'''
        write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
        txt132 = 'b_zgoubi.fai'
        write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
        txt132 = ' 1'
        write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
C Completes zgoubi.dat with the rest of zgoubi_StabLim-In.dat
 1      continue
          read(lunR,fmt='(a)',end=10) txt132
            if    (strcon(txt132,'''FAISTORE''',10,
     >                                        IS) ) then 
              read(lunR,fmt='(a)',end=62) txt132
              read(lunR,fmt='(a)',end=62) txt132
              goto 1
            elseif(strcon(txt132,'''FAISCNL''',9,
     >                                       IS) ) then 
              read(lunR,fmt='(a)',end=62) txt132
              goto 1
            elseif(strcon(txt132,'''TWISS''',7,
     >                                       IS) ) then 
              read(lunR,fmt='(a)',end=62) txt132
              goto 1
            elseif(strcon(txt132,'''CAVITE''',8,
     >                                       IS) ) then 
              read(lunR,fmt='(a)',end=62) txt132
              txt132 = txt132(debstr(txt132):finstr(txt132))
              if(txt132(1:2) .ne. txt0)
     >           txt132 = txt0//txt132(debstr(txt132):finstr(txt132))
              write(lunW,fmt='(a)')txt132(debstr(txt132):finstr(txt132))
              read(lunR,fmt='(a)',end=62) txt132
              write(lunW,fmt='(a)')txt132(debstr(txt132):finstr(txt132))   
              read(lunR,fmt='(a)',end=62) txt132
              write(lunW,fmt='(a)')txt132(debstr(txt132):finstr(txt132))
              goto 1
            endif
          if    (strcon(txt132,'''REBELOTE''',10,
     >                                        IS) ) then 
            write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))   
            read(lunR,fmt='(a)',end=62) txt132
            write(txt20,fmt='(I6)') nTurn
            txt132 = txt20//'  0.3  99'
            write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
            txt132 = '''END'''
            write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
            goto 10
          endif
          if    (strcon(txt132,'''END''',5,
     >                                     IS) ) then 
              txt132 = '''REBELOTE'''
              write(lunW,fmt='(a)')txt132(debstr(txt132):finstr(txt132))
              write(txt20,fmt='(i6)') nTurn
              txt132 = txt20//'  0.3  99'
              write(lunW,fmt='(a)')txt132(debstr(txt132):finstr(txt132))
              txt132 = '''END'''
              write(lunW,fmt='(a)')txt132(debstr(txt132):finstr(txt132))
              goto 10
          endif

          if(empty(txt132)) then
            write(lunW,fmt='(a)') '!'
          else
             write(lunW,fmt='(a)')
     >        txt132(debstr(txt132):finstr(txt132))
          endif
c          write(lunW,fmt='(a)') txt132   

        goto 1

 10     continue
        close(lunW)

c              write(*,*) ' done .... '
c                   stop

        call stabLim(HV,prec,jo,x(jo),xp(jo),z(jo),zp(jo),s(jo),d(jo),
     >                  xst,xpst,zst,zpst,sst,dst,duStrt,kex,ok)

        if(ok) then 
          jok = jok + 1
          write(6,*) 
          write(6,*) '---------  Found ',jok,' extreme stable orbits,',
     >      ' now #',jo,' :'
          write(6,fmt='(1p,4e14.6,e9.1,e12.4)') 
     >                          xst,xpst,zst,zpst,sst,dst
          write(6,*)
          write(icho,fmt='(a,i2,a,1p,e12.4,a,i2,a)') 
     >     ' Found ',jok,' extreme stable orbits, prec :',
     >      prec,' now #',jo,' :'
          write(icho,fmt='(1p,4e14.6,e9.1,e12.4)') 
     >                          xst,xpst,zst,zpst,sst,dst
          call flush2(icho,.false.)

          storb(1,jok) = xst
          storb(2,jok) = xpst
          storb(3,jok) = zst
          storb(4,jok) = zpst
          storb(5,jok) = sst
          storb(6,jok) = dst
          if(xst .gt. xma) xma = xst
          if(xpst .gt. xpma) xpma = xpst
          if(zst .gt. zma) zma = zst
          if(zpst .gt. zpma) zpma = zpst

        endif
        result = result .or. ok        

        write(*,*) ' jo / nTraj :    ', jo,' / ',nTraj
        if(jo.ge. nTraj) goto 60
        jo = jo+1

      goto 2

 60   continue
      close(lunR)

Create zgoubi_StabLim-Out.dat_H, or _Hz, or _V containing limit orbits
      namFil = 'zgoubi_StabLim-Out.dat'//'_'//HV
      open(unit=lunW,file=namFil)
      open(unit=lunR,file='zgoubi_StabLim-In.dat')
C Read/write till "KOBJ"
        txt132 = ' '
        dowhile
     >  (txt132(2:6) .ne. 'OBJET' .and. TXT132(2:8) .ne. 'MCOBJET')
          read(lunR,fmt='(a)') txt132
          txt132 = txt132(debstr(txt132):finstr(txt132))
          if(empty(txt132)) then
            write(lunW,fmt='(a)') '!'
          else
             write(lunW,fmt='(a)')
     >        txt132(debstr(txt132):finstr(txt132))
          endif
c          write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
        enddo
        do i=1,2
          read(lunR,fmt='(a)') txt132
          if(empty(txt132)) then
            write(lunW,fmt='(a)') '!'
          else
             write(lunW,fmt='(a)')
     >        txt132(debstr(txt132):finstr(txt132))
          endif
c          write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
        enddo
C Read/write "IMAX IMAXT"
        read(lunR,fmt='(a)') txt132
        if(kprx.eq.1) then
          kprx1 = 1
        else
          kprx1 = 0
        endif
        write(lunW,fmt='(i3,a)') kprx1+jok*njj,'  1'
        if(kprx.eq.1) write(lunW,fmt='(a)') 
     >    ' .0001  .0001  .0001  .0001  .0001  1.  ''o'''
C Write stab Lim coordinates
        do j=1,jok

          if    (njj .eq. 1) then
            dxx = 0.d0 ;       dzz = 0.d0
            xxx = storb(1,j) ; zzz = storb(3,j) 
          elseif(njj .ge. 2) then
            dxx = (storb(1,j) - xi(j))/float(njj-1)
            dzz = (storb(3,j) - zi(j))/float(njj-1)
            xxx = xi(j) -dxx
            zzz = zi(j) -dzz
          else
            write(*,*) 
            write(*,*) 'Pgm searStabLim.  Found njj = ',njj,'  !!'
            stop 'Such njj value cannot be !'
          endif

c         write(*,*) ' initial : ',xi(j),zi(j)
c         write(*,*) ' xxx, dxx ',storb(1,j),xi(j),xxx-dxx, dxx
c         write(*,*) ' zzz, dzz ',storb(3,j),zi(j),zzz-dzz, dzz
c                      read(*,*)

         read(lunR,fmt='(a)') txt132
          do jj = 1, njj
            if    (HV .eq. 'H' .or. HV .eq. 'Hz') then
C              xxx = storb(1,j)/float(njj) * float(jj)
              xxx = xxx + dxx
              zzz = storb(3,j)
            elseif(HV .eq. 'V') then
              xxx = storb(1,j)
C              zzz = storb(3,j)/float(njj) * float(jj)
              zzz = zzz + dzz
            endif
            write(lunW,fmt='(1p,4e16.8,e9.1,e16.8,4a,i4)') 
     >      xxx, storb(2,j), zzz, storb(4,j),
     >      storb(5,j), storb(6,j),' ','''',let(j),''' ',j
          enddo
        enddo
        do j=1,nTraj-jok
          read(lunR,fmt='(a)') txt132
        enddo
C Complete OBJET with the line of 1's
 33     read(lunR,fmt='(a)',end=62) txt132
        txt132 = txt132(debstr(txt132):finstr(txt132))
        if(txt132(1:2) .eq. '1 ') then
          goto 33
        else
          backspace(lunR)
        endif
        ii = 0
 34     write(lunW,*) ' 1 1 1 1 1 1 1 1 1 1 ' 
        ii = ii + 10
        if(kprx1+jok*njj .gt. ii) goto 34
        txt132 = '''FAISTORE'''
        write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
        txt132 = 'b_zgoubi.fai'
        write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
        txt132 = ' 1'
        write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
C Completes zgoubi_StabLim-Out.dat_HV with the rest of zgoubi_StabLim-In.dat
 11     continue
          read(lunR,fmt='(a)',end=62) txt132

          if    (strcon(txt132,'''FAISTORE''',10,
     >                                        IS) ) then 
            read(lunR,fmt='(a)',end=62) txt132
            read(lunR,fmt='(a)',end=62) txt132
            goto 11
          elseif(strcon(txt132,'''FAISCNL''',9,
     >                                       IS) ) then 
            read(lunR,fmt='(a)',end=62) txt132
            goto 11
          elseif(strcon(txt132,'''TWISS''',7,
     >                                       IS) ) then 
            read(lunR,fmt='(a)',end=62) txt132
            goto 11
          elseif(strcon(txt132,'''CAVITE''',8,
     >                                       IS) ) then 
            read(lunR,fmt='(a)',end=62) txt132
            txt132 = txt132(debstr(txt132):finstr(txt132))
            if(txt132(1:2) .ne. txt0) 
     >         txt132 = txt0//txt132(debstr(txt132):finstr(txt132))
            write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
            read(lunR,fmt='(a)',end=62) txt132
            write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))   
            read(lunR,fmt='(a)',end=62) txt132
            write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))   
            goto 11
          elseif(strcon(txt132,'''REBELOTE''',10,
     >                                        IS) ) then 
            read(lunR,fmt='(a)',end=62) txt132
            goto 11
          elseif    (strcon(txt132,'''END''',5,
     >                                     IS) ) then 
            goto 63
          endif

          if(empty(txt132)) then
            write(lunW,fmt='(a)') '!'
          else
             write(lunW,fmt='(a)')
     >        txt132(debstr(txt132):finstr(txt132))
          endif
c          write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))   

        goto 11

 63   CONTINUE
      txt132 = '''REBELOTE'''
      write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
      write(txt20,fmt='(I6)') nTurn
      txt132 = txt20//'  0.3  99'
      write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
      txt132 = '''END'''
      write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132)) 
 62   CONTINUE
      close(lunR)
      close(lunW)

      write(*,*) ' '
      write(*,*) '--------------'
      write(*,fmt='(3a,i2,a,i2,a)') ' New file ',namFil,
     >' contains the ',jok,
     >' traj. below (over ',nTraj,' trajectories launched) :'
      do j=1,jok
        write(*,fmt='(1p,4e14.6,e9.1,e12.4,4a,i4)') 
     >  storb(1,j), storb(2,j), storb(3,j), storb(4,j),
     >  storb(5,j), storb(6,j),' ','''',let(j),''' ',j
      enddo
      write(*,*) ' '
      if(first) then
        call system('cat searchStabLim.out >> searchStabLim.out_old')
        first = .false.
      endif

      open(unit=33,file='searchStabLim.temp')
      WRITE(33,*) 
     >   '% Stab. limit coordinates 1-6, tag, traj#, K, xi'
      do j=1,jok
        write(33,fmt='(1p,4e14.6, e9.1, e12.4, 4a, i4, a,a2)') 
     >  storb(1,j), storb(2,j), storb(3,j), storb(4,j),
     >  storb(5,j), storb(6,j),' ','''',let(j),''' ',j ,' ',HV
      enddo
      write(33,*) ' '
      close(33)
      call system('cat searchStabLim.temp >> searchStabLim.out')
      call system('rm searchStabLim.temp')

          open(unit=34,file='tempHV.dum')
          write(34,*) HV
          close(34)
      if    (HV .eq. 'V') then
          call system('cp -f zgoubi_StabLim-Out.dat_V zgoubi.dat')
      elseif(HV .eq. 'H') then
          call system('cp -f zgoubi_StabLim-Out.dat_H zgoubi.dat')
      elseif(HV .eq. 'Hz') then
          call system('cp -f zgoubi_StabLim-Out.dat_Hz zgoubi.dat')
      endif

      if(kTrk.eq.1) then
        write(*,*)  '  '
        write(*,*)  '///////////////////////////////////////'
        write(*,*)  ' Now tracking for tunes '
        write(*,*)  '///////////////////////////////////////'
        write(*,*)  '  '
        call system('~/zgoubi/source/zgoubi')
        call system('~/zgoubi/toolbox/tunesFromFai/tunesFromFai')  
        if    (HV .eq. 'V') then
          call system('mv -f b_zgoubi.fai b_zgoubi.fai_V')
          call system('mv -f tunesFromFai.out tunesFromFai.out_V')
        elseif(HV .eq. 'H') then
          call system('mv -f b_zgoubi.fai b_zgoubi.fai_H')
          call system('mv -f tunesFromFai.out tunesFromFai.out_H')
        elseif(HV .eq. 'Hz') then
          call system('mv -f b_zgoubi.fai b_zgoubi.fai_Hz')
          call system('mv -f tunesFromFai.out tunesFromFai.out_Hz')
        endif

      endif

      goto 111   
C----- Loop on HV

c      write(*,*) ' SearchStabLim :  end of loop 111... '
      goto 99

 997  continue
      write(*,*) ' SearchStabLim :  end of job on 997... '
      goto 99

 999  continue
      write(*,*) ' *** SearchStabLim : error upon reading lunData '
      goto 99

 998  continue
      write(*,*) ' *** SearchStabLim : end of file lunData reached'
      goto 99

 99   continue
      write(*,*) ' JOB TERMINATED '
      close(lunData)
      stop
      end
      FUNCTION STRCON(STR,STRIN,NCHAR,
     >                                IS)
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
      subroutine stabLim(HV,prec,jo,xco,xpco,zco,zpco,s,d,
     >                      xst,xpst,zst,zpst,sst,dst,duStrt,kex,ok)
      character*(*) HV
      logical lost, ok, first
      parameter (lunR=13,lunW=14)
      character txt132*132, let*1
      data first / .true. /
      data let / 'i' /
      logical du2
      integer debstr, finstr
      logical  empty

      ok = .false.
      lost = .false.
      du2 = .false.

C du stands for dx or dz (cm)
C      du = 5.   !!cm 
C      du = 50.   !!cm 
      du = duStrt   !!cm 
          if(HV .eq. 'V') then
            zco = zco + du
            zst = zco 
            xst = xco
          elseif(HV .eq. 'H' .or. HV .eq. 'Hz') then
            xco = xco + du
            xst = xco 
            zst = zco
          endif
      xpst = xpco
      zpst = zpco
      sst = s
      dst = d

 2    continue

C----------------------- Rebuild zgoubi.dat with new traj. -------------
        call system('cp -f zgoubi.dat searchStabLim.temp2')
        close(lunR)
        close(lunW)
        open(unit=lunR,file='searchStabLim.temp2')
        open(unit=lunW,file='zgoubi.dat')
C Read/write till "KOBJ"
        txt132 = ' '
        dowhile
     >  (txt132(2:6) .ne. 'OBJET' .and. TXT132(2:8) .ne. 'MCOBJET')
          read(lunR,fmt='(a)') txt132
          txt132 = txt132(debstr(txt132):finstr(txt132))
          if(empty(txt132)) then
            write(lunW,fmt='(a)') '!'
          else
             write(lunW,fmt='(a)')
     >        txt132(debstr(txt132):finstr(txt132))
          endif
c          write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
        enddo
          do i=1,2
            read(lunR,fmt='(a)') txt132
          if(empty(txt132)) then
            write(lunW,fmt='(a)') '!'
          else
             write(lunW,fmt='(a)')
     >        txt132(debstr(txt132):finstr(txt132))
          endif
c            write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))
          enddo
C Read/write "IMAX IMAXT"
          read(lunR,fmt='(a)') txt132
          write(lunW,fmt='(a)') '1  1'
C Write best co coordinates
          read(lunR,fmt='(a)') txt132
          write(lunW,fmt='(1p,4e14.6,e9.1,e12.4,4a,i4)') 
     >      xst,xpst,zst,zpst,sst,dst,' ','''',let,''' ',jo
          write(*,fmt='(1p,a,4e14.6,e9.1,e12.4)') 
     >    'Initial object :', xco,xpco,zco,zpco,sst,dst
          write(*,fmt='(1p,a,4e14.6,e9.1,e12.4,4a,i4)') 
     >    'New object :', xst,xpst,zst,zpst,sst,dst,' ','''','A','''',jo
C Completes zgoubi.dat with the rest of searchStabLim.temp2
 19       continue
          read(lunR,fmt='(a)',end=10) txt132
          if(empty(txt132)) then
            write(lunW,fmt='(a)') '!'
          else
             write(lunW,fmt='(a)')
     >        txt132(debstr(txt132):finstr(txt132))
          endif
c          write(lunW,fmt='(a)') txt132(debstr(txt132):finstr(txt132))   
        goto 19 
 10     continue
        close(lunR)
        close(lunW)
C-----------------------------------------------------------------------

C Run zgoubi with single traj
        call system('~/zgoubi/source/zgoubi')
    
C Get initial Traj coord from zgoubi.res
        call getRes(
     >              lost)

C      write(*,*) ' x,xp,z,zp,s,d,lost : ',xst,xpst,zst,zpst,sst,dst,lost

        if(lost) then
          write(6,*) 
          write(6,*) '-----  Particle # ',jo,' lost, back to previous '
          write(6,*) '         du -> -du =',-du
          write(6,*)

          if(HV .eq. 'V') then
            ylost = zco
            zco = zco - du
            zst = zco 
            xst = xco
          elseif(HV .eq. 'H' .or. HV .eq. 'Hz') then
            ylost = xco
            xco = xco - du
            xst = xco 
            zst = zco
          endif

          du = du/2.
          du2 = .true.
          if (abs(du).le.prec) goto 98
          lost = .false.
          call system('mv -f b_zgoubi.fai b_zgoubi-temp.fai')

        endif      

        write(6,*) '-----  Particle # ',jo,', increase initial ',
     >                                         HV,' position...'
        if(du2) du = du/2.
        du2 = .false.
        write(6,*) '  du = ', du
          if(HV .eq. 'V') then
            zco = zco + du
            if(zco.eq.ylost) goto 98
            zst = zco 
            xst = xco
          elseif(HV .eq. 'H' .or. HV .eq. 'Hz') then
            xco = xco + du
            if(xco.eq.ylost) goto 98
            xst = xco 
            zst = zco
          endif

        iter = iter+1
c        if(iter.le.10) then
           goto 2
c        else
c           ok = .false.
c           goto 99
c        endif

 98        continue
            ok = .true.
            lost = .true.
            call system('mv -f b_zgoubi-temp.fai b_zgoubi.fai')
            goto 99

 99   continue
      write(*,*) ' x,xp,z,zp,du,lost : ',xst,xpst,zst,zpst,du,lost
      call system('rm searchStabLim.temp2')
      return
      end
      subroutine getRes(
     >                  lost)
C Get initial Traj coord, and average orbit, from zgoubi.res
      character let*1
      logical lost
C Read pick-ups from last pass in zgoubi.res, and cumulates in readPU.out
C This software assumes use of REBELOTE with writes inhibited - so that PU readouts 
C are written in zgoubi.res at first and last pass only. 

      character txt132*132

      logical strcon
      integer debstr, finstr
      data lunR / 15 /

      open(unit=lunR,file='zgoubi.res') 

 12   continue
      read(lunR,fmt='(a)',end=18,err=19) txt132

      if(strcon(txt132,'lost',4,
     >                          IS)) then
          goto 11
      else
        goto 12
      endif

 11   continue
        write(*,*) '****  Found  part. lost, with x, z =',x, z
        write(6,*)  ' ' 
        lost = .true.
      goto 99

 18   continue
      lost = .false.
      write(*,*) ' sbr getRes * End of read, no loss *   '
      goto 99

 19   continue
      lost = .false.
      write(*,*) ' sbr getRes * End of readPU upon read error *   '

 99   close(lunR)
      return
      end
      subroutine gotoEnd(lunOut,
     >                          nTrajRun)
      character*20 txt
      logical strcon

        write(6,*) '  Now go to End of readPU.out'

      nTrajRun = 0
C On cherche 'Closed'
 10   read(lunOut,fmt='(A)',end=98,err=99) txt
C      write(6,*) ' txt : ', txt
      if(strcon(txt,'% Closed',8,
     >                         IS)) then 
        nTrajRun = nTrajRun + 1 
C        write(6,*) '   gotoEnd,  nTrajRun = ', nTrajRun
      endif
      goto 10

      return

 98   write(6,*) '  * End of readPU/gotoEnd upon EOFile lunOut *    '
      return
 99   write(6,*) '  * End of readPU/gotoEnd upon read error in lunOut *'
      return
      end
      SUBROUTINE FLUSH2(IUNIT,BINARY)
      LOGICAL BINARY
      CHARACTER*80 TXT80
      BACKSPACE(IUNIT)
      IF(.NOT.BINARY) THEN
        READ(IUNIT,FMT='(A80)') TXT80
      ELSE
        READ(IUNIT) TXT80
      ENDIF
      RETURN
      END
      FUNCTION EMPTY(STR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EMPTY
      CHARACTER*(*) STR
C     -----------------------------------------------------
C     .TRUE. if STR is either empty or contains only blanks
C     -----------------------------------------------------

      INTEGER FINSTR
      EMPTY = FINSTR(STR) .EQ. 0
      RETURN
      END
