C Starting from initial coordinates near closed orbits, as read from input file 
C zgoubi.dat, search stability limits upon scan of x,z coordinates. 
C zgoubi.dat is presumed to contain stable trajectories, which will be pushed to 
C extreme amplitude
      parameter (lunR=11,lunW=12,lunIn=15,lunSto=18)
      character txt132*132, txt6*6, let*1, namFil*30
      parameter (nTrajmx=99)
      dimension x(nTrajmx),xp(nTrajmx),z(nTrajmx),
     >                   zp(nTrajmx),s(nTrajmx),d(nTrajmx),let(nTrajmx)
      dimension storb(7,nTrajmx)
      logical ok, result, strcon, first
      logical okFAI, okREB
      data okFAI / .true. /
C H=pure H, Hz=H+small z, V=vertical 
      data lunData / 7 /
C      prec : cm
      data prec, nTurn / 2.,  29 /
      data ok, first, okREB / .false., .true., .false. /

C---------- Input data ------------------
      open(unit=lunData,file='searchDA.data')
      read(lunData,*,err=66,end=66) precIn,dxIn,dzIn,xMaxIn,nTrnIn
      close(lunData)
      prec = precIn
      dx = dxIn
      dz = dzIn
      xMax = xMaxIn
      nTurn = nTrnIn
      goto 67
 66   continue
      write(6,*) '  Could not read from searchDA.data. Default values, '      
 67   continue
      write(6,*) '  prec, dz, nTurn  : ',prec,dz,nTurn
C----------------------------------------
      call system('cat searchDA.out >>  searchDA.out_old ; 
     >             rm -f searchDA.out')
      open(unit=lunSto,file='searchDA.out')
          open(unit=35,file='tempKXi.dum')
          read(35,*) xK,xiDeg
          close(35)
          write(lunSto,*) '%  xK,xiDeg : ',xK,xiDeg

      open(unit=lunData,file='searchDA.data')

      result = .false.

      call system('rm -f zgoubi.dat')

C zgoubi_searchDA-In.dat is supposed to contain stable trajectories for various momentum values
C - these could be output from prior searchCO procedure - 
C these will be starting point for search of DA
C
      open(unit=lunR,file='zgoubi_searchDA-In.dat')

      jo = 1
      jok = 0
 2    continue
 
        open(unit=lunW,file='zgoubi.dat')
        rewind(lunW)
        rewind(lunR)

C Read till "KOBJ"
        do i=1,4
          read(lunR,fmt='(a)') txt132
          write(lunW,*) txt132
        enddo
C Read till "IMAX IMAXT"
          read(lunR,*) nTraj, imaxt
          if(nTraj.gt.nTrajmx) stop ' Too many trajectories...'
          txt132 = ' 1  1'
          write(lunW,*) txt132
C Read all initial traj from zgoubi_searchDA-In.dat, supposed to be stable
        do i=1,nTraj
          read(lunR,*) x(i),xp(i),z(i),zp(i),s(i),d(i),let(i)
        enddo
C Retains only one initial traj at a time
C and sets z to nul
        z(jo) = 0.
        zp(jo) = 0.
        write(lunW,fmt='(1p,4e14.6,e9.1,e12.4,4a,i4)') 
     >    x(jo),xp(jo),z(jo),zp(jo),s(jo),d(jo),' ','''',let(jo),'''',jo
        write(6,*) 
        write(6,*) ' ---------    New trajectory :'
        write(6,fmt='(1p,4e14.6,e9.1,e12.4,4a,i4)') 
     >    x(jo),xp(jo),z(jo),zp(jo),s(jo),d(jo),' ','''',let(jo),'''',jo
        write(6,*)
C Completes zgoubi.dat with the rest of zgoubi_searchDA-In.dat
 1      continue
          read(lunR,fmt='(a)',end=10) txt132
C           skip storage, so to save on CPU time :

          if(.not. okFAI) then
            if    (strcon(txt132,'FAISTORE',8,
     >                                        IS) ) then 
              read(lunR,fmt='(a)',end=62) txt132
              read(lunR,fmt='(a)',end=62) txt132
              read(lunR,fmt='(a)',end=62) txt132
            elseif(strcon(txt132,'FAISCNL',7,
     >                                       IS) ) then 
              read(lunR,fmt='(a)',end=62) txt132
              read(lunR,fmt='(a)',end=62) txt132
            endif
          endif
          if    (strcon(txt132,'REBELOTE',8,
     >                                        IS) ) then 
            okREB = .true.
            write(lunW,*) txt132   
            read(lunR,fmt='(a)',end=62) txt132
            write(txt6,*) nTurn
            txt132 = txt6//'  0.2  99'
            write(lunW,*) txt132
            txt132 = '''END'''
            write(lunW,*) txt132
            goto 10
          endif
          if    (strcon(txt132,'''END''',5,
     >                                     IS) ) then 
            if(.not. okREB) then
              txt132 = '''REBELOTE'''
              write(lunW,*) txt132
              write(txt6,*) nTurn
              txt132 = txt6//'  0.2  99'
              write(lunW,*) txt132
              txt132 = '''END'''
              write(lunW,*) txt132
              goto 10
            endif
          endif

          write(lunW,*) txt132   

        goto 1

 10     continue
        close(lunW)

        call DA(
     >    lunSto,jo,dx,dz,xMax,prec,x(jo),xp(jo),z(jo),zp(jo),
     >    s(jo),d(jo),xst,xpst,zst,zpst,sst,dst,kex,ok)

        if(ok) then 

          jok = jok + 1
          write(6,*) 
          write(6,*) '---------  Found ',jok,' extreme stable orbits,',
     >      ' now #',jo,' :'
          write(6,fmt='(1p,4e14.6,e9.1,e12.4)') 
     >                          xst,xpst,zst,zpst,sst,dst
          write(6,*)

          storb(1,jok) = xstP
          storb(7,jok) = xstM
          storb(2,jok) = xpst
          storb(3,jok) = zst
          storb(4,jok) = zpst
          storb(5,jok) = sst
          storb(6,jok) = dst
        endif
        result = result .or. ok        

        write(*,*) ' jo / nTraj :    ', jo,' / ',nTraj
        if(jo.ge. nTraj) goto 60
        jo = jo+1

      goto 2

 60   continue
 62   continue
      close(lunR)
      close(lunW)

      goto 99

 999  continue
      write(*,*) ' *** Error upon reading lunData '
      goto 99

 998  continue
      write(*,*) ' *** End of file lunData reached'
      goto 99

 99   continue
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
      subroutine DA(lun,jo,dx1,dz1,xMax,prec,xco,xpco,zco,zpco,
     >                 s,d,xst,xpst,zst,zpst,sst,dst,kex,ok)
      logical lost, ok, first
      parameter (lunR=13,lunW=14)
      character txt132*132, let*1
      data first / .true. /
      data let / 'i' /
      parameter (zero=0.)
      character*2 case

      ok = .false.
      lost = .false.

      xst = xco - xMax
      xpst = xpco
      zst = 0.
      zpst = 0.
      sst = s
      dst = d

      dx = dx1
      dz = dz1
      iz = 1
      ix = 1

 2    continue

C----------------------- Rebuild zgoubi.dat with new traj. -------------
        call system('cp -f zgoubi.dat searchDA.temp2')
        close(lunR)
        close(lunW)
        open(unit=lunR,file='searchDA.temp2')
        open(unit=lunW,file='zgoubi.dat')
C Read/write till "KOBJ"
          do i=1,4
            read(lunR,fmt='(a)') txt132
            write(lunW,fmt='(a)') txt132
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
C Completes zgoubi.dat with the rest of searchDA.temp2
 19       continue
          read(lunR,fmt='(a)',end=10) txt132
          write(lunW,fmt='(a)') txt132   
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

        if(lost) then
          if(ix.eq.1) then
            if(abs(dx).le.prec) then
C No more hope !           
              ok = .true.
              xlostP = xst
              xlostM = xst
              zlost = zst
                write(lun,fmt='(1p,2e14.6,2i4,a)')   
     >                xlostM,zlost,iz,ix, ' xlost_min, zlost, iz, ix'
                write(lun,fmt='(1p,2e14.6,2i4,a)') 
     >                xlostP,zlost,iz,ix, ' xlost_max,   "  , " , " '
                CALL FLUSH2(lun,.FALSE.)
                goto 99
            else
                ix = 0
                dx = -dx
            endif
          else
            if(case.eq.'xP') then
C              x est croissant a partir xco 
              if(dx.le.prec) then
C                => on inverse pour l'autre coté   
                xlostP = xst
                xst = xco - (xlostP-xco)
                dx = - (xlostP-xco)/10.
                case = 'xM'
              else
                xst = xst -dx
                dx = dx/2.
              endif
            else  ! case = 'xM'
C             les deux sens en x on été scannés ; on incremente z
              if(-dx.le.prec) then
                xlostM = xst
                zlost = zst
                write(lun,fmt='(1p,2e14.6,2i4,a)')   
     >                xlostM,zlost,iz,ix, ' xlost_min, zlost, iz, ix'
                write(lun,fmt='(1p,2e14.6,2i4,a)') 
     >                xlostP,zlost,iz,ix, ' xlost_max,   "  , " , " '
                CALL FLUSH2(lun,.FALSE.)
                dx = (xlostP-xco)/10.
                xst = xlostP
                zst = zst + dz
                iz = iz + 1
                ix = 0
                case = 'xP'
              else
                xst = xst - dx
                dx = dx / 2.
              endif
            endif
            write(*,*)' xst, zst, dx, dz, iz, ix',xst,zst,dx,dz,iz,ix
          endif
        endif      

        ix = ix+1
        xst = xst + dx
        call system('mv -f b_zgoubi.fai b_zgoubi-temp.fai')

        goto 2

 99   continue
      write(*,*) ' x,xp,z,zp,dx,dz,lost : ',xst,xpst,zst,zpst,dx,dz,lost
      call system('rm searchDA.temp2')
      return
      end
      subroutine DA_Old(jo,dx,dz,xco,xpco,zco,zpco,s,d,
     >                      xst,xpst,zst,zpst,sst,dst,kex,ok)
      logical lost, ok, first
      parameter (lunR=13,lunW=14)
      character txt132*132, let*1
      data first / .true. /
      data let / 'i' /
      parameter (zero=0.)

      ok = .false.
      lost = .false.

      xst = xco 
      xpst = xpco
      zst = zco
      zpst = zpco
      sst = s
      dst = d

      iz = 1
      ix = 1

 2    continue

C----------------------- Rebuild zgoubi.dat with new traj. -------------
        call system('cp -f zgoubi.dat searchDA.temp2')
        close(lunR)
        close(lunW)
        open(unit=lunR,file='searchDA.temp2')
        open(unit=lunW,file='zgoubi.dat')
C Read/write till "KOBJ"
          do i=1,4
            read(lunR,fmt='(a)') txt132
            write(lunW,fmt='(a)') txt132
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
C Completes zgoubi.dat with the rest of searchDA.temp2
 19       continue
          read(lunR,fmt='(a)',end=10) txt132
          write(lunW,fmt='(a)') txt132   
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

        if(lost) then
          if(xst .eq. xco) then
C No more hope !           
            ok = .true.
            xlostP = xst
            xlostM = xst
            zlost = zst
            write(88,fmt='(1p,2e14.6,i4)')   xlostM,zlost,iz
            write(88,fmt='(1p,2e14.6,i4)')   xlostP,zlost,iz
            call system('mv -f b_zgoubi-temp.fai b_zgoubi.fai')
            goto 99
          else
            if(dx.gt.0) then
C              x est croissant a partir xco => on inverse pour l'autre coté               
              xlostP = xst
              dx = -dx
              xst = xco + dx
            else
C             les deux sens en x on été scannés ; on incremente z
              xlostM = xst
              zlost = zst

              write(88,fmt='(1p,2e14.6,i4)')   xlostM,zlost,iz
              write(88,fmt='(1p,2e14.6,i4)')   xlostP,zlost,iz

              xst = xco
              dx = -dx
              zst = zst + dz
              iz = iz + 1
            endif
            write(*,*) ' xst, zst, dx, dz, iz',  xst, zst, dx, dz, iz
C            pause
          endif
        else
          xst = xst + dx
          call system('mv -f b_zgoubi.fai b_zgoubi-temp.fai')
        endif      

        goto 2

 99   continue
      write(*,*) ' x,xp,z,zp,dx,dz,lost : ',xst,xpst,zst,zpst,dx,dz,lost
      call system('rm searchDA.temp2')
      return
      end
      subroutine getRes(
     >                  lost)
C Get initial Traj coord, and average orbit, from zgoubi.res
      character let*(*)
      logical lost
C Read pick-ups from last pass in zgoubi.res, and cumulates in readPU.out
C This software assumes use of REBELOTE with writes inhibited - so that PU readouts 
C are written in zgoubi.res at first and last pass only. 

      character txt132*132

      logical strcon
      integer debstr, finstr
      data lunR / 15 /

      open(unit=lunR,name='zgoubi.res') 

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
      subroutine gotoEnd(lun)
      character*132 txt
      iter = 0
 10   read(lunOut,fmt='(A)',end=98,err=98) txt
      iter = iter+1
      goto 10
 98   write(6,*) '  * End of searchDA.out reached afetr ',iter,' lines'
      write(lun,*) '   '
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