C Starting from initial coordinates near closed orbits, at various momenta,  
C as read from input file zgoubi.dat, search z-stability limit 
C for a series of x=xco,xco+k*dx. 
C zgoubi.dat is presumed to contain stable trajectories, which will be pushed to 
C extreme amplitude
      parameter (lunR=11,lunW=12,lunIn=15,lunSto=18,icho=10)
      character txt132*132, txt6*6, let*1, namFil*30
      parameter (nTrajmx=99)
      dimension x(nTrajmx),xp(nTrajmx),z(nTrajmx),
     >                   zp(nTrajmx),s(nTrajmx),d(nTrajmx),let(nTrajmx)
C      dimension dxminu(nTrajmx), dxplus(nTrajmx), duStrt(nTrajmx)
      dimension xminu(nTrajmx), xplus(nTrajmx), duStrt(nTrajmx)
      parameter (ixMax=400)
      dimension storb(7,ixMax)
      logical ok, strcon, first
      character*2 HV, HVIn
      logical   chs

      INTEGER DEBSTR,FINSTR

C H=pure H, Hz=H+small z, V=vertical 
      data HV / 'V' /
      data lunData / 7 /
C                          prec : cm
      data zDA, prec, nTurn / 0., 2., 29 /  
      data ok, first/ .false., .true. /
C------------------------------
C These data may be changed
      data klostM / 2 /
C       duStrt : cm
!      data duStrt / nTrajmx*50. /  ! Good for NuFact storage ring
      data duStrt / nTrajmx*2. /  ! Good for superB
C------------------------------
      open(unit=icho,file='searchDA.echo')

C---------- Input data ------------------
      open(unit=lunData,file='searchDA.data')
C      read(lunData,*,err=66,end=66) NdxIn,preczIn,nTrnIn,duStrtIn
      read(lunData,*,err=66,end=66) NdxIn,preczIn,nTrnIn
      close(lunData)
      prec = preczIn
      Ndx = NdxIn
      nTurn = nTrnIn
      nTr12 = nTurn  ! /2 + 1
C      duStrt = duStrtIn
      goto 67
 66   continue
      write(6,*) '  Could not read from searchDA.data'
      stop  '  Check content of  searchDA.data '
 67   continue
      write(6,*) '  #dx, z-prec, nTurn  : ',Ndx, prec,nTurn
      write(icho,fmt='(2a,i4,1x,1p,e12.4,1x,i4)') 'Read searchDA.data ;'
     > ,' #dx, z-prec, nTurn  : ',Ndx, prec,nTurn
      call flush2(icho,.false.)

C----------------------------------------
C First search min x-, max x+  so to determine dx-, dx+, and max z to get starting z-value
      call system('cp zgoubi_StabLim-In.dat zgoubi_StabLim-In.dat-save')
      call system('cp zgoubi_searchDA-In.dat zgoubi_StabLim-In.dat')

! Stage 2 :  search lower stable limit x-
!   2.a create searchStabLim.data for H lower limit x- (including subliminal z)
      write(icho,fmt='(a)') ' Now doing StabLim H-'
      call flush2(icho,.false.)
      open(unit=lunw,file='searchStabLim.data')
      write(lunw,fmt='(a,1p,1(e14.6,1x),a,i5,a)')
     >'H ',prec,' -1. ',nTr12,' 0 1 0'
C      write(lunw,fmt='(a,1p,2(e14.6,1x),i5,a)')
C     >'H ',prec,-abs(duStrt),nTr12,' 0 1 0'    ! cyclotron Luciano
      write(lunw,*) 'H[z]/V, prec(cm), strtng dx or dz(cm),'
     >,' #turns, final trckng y/n=1/0, #lips to trck, prxial lips,'
     >,' (y/n=1/0)'

      close(lunw)
!   2.b search x-. This will also create zgoubi_StabLim-Out.dat_H
      call system(
     >'~/zgoubi/struct/tools/searchStabLim/searchStabLim')
!   2.c Read x-
      open(unit=lunw,file='zgoubi_StabLim-Out.dat_H')
      do i = 1, 5
        read(lunw,fmt='(a132)') txt132
      enddo
      read(txt132,*) ntraj
      do jo = 1, ntraj
        read(lunw,*) temp
        xminu(jo) = temp 
c        dxminu(jo) = abs(xminu) / float(Ndx)
      enddo
      close(lunw)
      call system('mv zgoubi_StabLim-Out.dat_H 
     >      zgoubi_StabLim-Out.dat_Hm')

! Stage 3 : search upper stable limit x+
!   3.a create searchStabLim.data for H upper limit x+ (including subliminal z)
      write(icho,fmt='(a)') ' Now doing StabLim H+'
      call flush2(icho,.false.)
      open(unit=lunw,file='searchStabLim.data')
      write(lunw,fmt='(a,1p,1(e14.6,1x),a,i5,a)')
     >'H ',prec,' +1. ',nTr12,' 0 1 0'
C      write(lunw,fmt='(a,1p,2(e14.6,1x),i5,a)')
C     >'H ',prec, abs(duStrt),nTr12,' 0 1 0'    ! cyclotron Luciano
      write(lunw,*) 'H[z]/V, prec(cm), strtng dx or dz(cm),'
     >,' #turns, final trckng y/n=1/0, #lips to trck, prxial lips,'
     >,' (y/n=1/0)'
      close(lunw)
!   3.b search x+. This will also create zgoubi_StabLim-Out.dat_H
      call system(
     >'~/zgoubi/struct/tools/searchStabLim/searchStabLim')
!   3.c Read x+
      open(unit=lunw,file='zgoubi_StabLim-Out.dat_H')
      do i = 1, 5
        read(lunw,fmt='(a132)') txt132
      enddo
      read(txt132,*) ntraj
      do jo = 1, ntraj
        read(lunw,*) temp
        xplus(jo) = temp 
C        dxplus(jo) = abs(xplus) / float(Ndx)
      enddo
      close(lunw)
      call system('mv zgoubi_StabLim-Out.dat_H 
     >      zgoubi_StabLim-Out.dat_Hp')

! Stage 4 :  search  z-limit at x=0
!   4.a create searchStabLim.data for V  limit z0
      write(icho,fmt='(a)') ' Now doing StabLim V at all x0'
      call flush2(icho,.false.)
      open(unit=lunw,file='searchStabLim.data')
      write(lunw,fmt='(a,1p,1(e14.6,1x),a,i5,a)')
     >'V ',prec,' +1. ',nTr12,' 0 1 0'
C      write(lunw,fmt='(a,1p,2(e14.6,1x),i5,a)')
C     >'V ',prec, abs(duStrt),nTr12,' 0 1 0'    ! cyclotron Luciano
      write(lunw,*) 'H[z]/V, prec(cm), strtng dx or dz(cm),'
     >,' #turns, final trckng y/n=1/0, #lips to trck, prxial lips,'
     >,' (y/n=1/0)'
      close(lunw)
!   4.b search z0. This will also create zgoubi_StabLim-Out.dat_V
      call system(
     >'~/zgoubi/struct/tools/searchStabLim/searchStabLim')
!   4.c Read z0
      open(unit=lunw,file='zgoubi_StabLim-Out.dat_V')
      do i = 1, 5
        read(lunw,fmt='(a132)') txt132
      enddo
      read(txt132,*) ntraj
      do jo = 1, ntraj
        read(lunw,*) bid, bid, z0
        if(z0 .ge. prec) then
          duStrt(jo) = 1.9*z0
        else
          duStrt(jo) = abs(xplus(jo))
          if(abs(xminu(jo)).gt.abs(xplus(jo))) duStrt(jo)=abs(xminu(jo))
        endif
      enddo
      close(lunw)
!!!      call system('rm zgoubi_StabLim-Out.dat_V')
      call system('rm zgoubi_StabLim-In.dat')

C----------------------------------------------
      write(icho,fmt='(a)') ' Now start searchDA'
      call flush2(icho,.false.)

      call system('cat searchDA.out >>  searchDA.out_old ; 
     >             rm -f searchDA.out')
      open(unit=lunSto,file='searchDA.out')

      do jo = 1, ntraj
        write(lunSto,fmt='(1p,a,3(e12.4,1x),1x,i3)')
     >  ' dx-, dx+, du_Z, itraj : ',xminu(jo),xplus(jo),duStrt(jo),jo
      enddo
      write(lunSto,*) '#  prec, Ndx, nTrnIn : ',prec,Ndx,nTrnIn
      write(lunSto,*) '# '
      write(lunSto,*) '# xst-xst0, zDA, dpp, xst0, jo, klost'
      write(lunSto,*) '# '
      call flush2(lunSto,.false.)

      call system('rm -f zgoubi.dat')

C zgoubi_searchDA-In.dat is supposed to contain stable trajectories for various momentum values
C - these could be output from prior searchCO procedure - 
C these will be starting point for search of DA
      
      open(unit=lunR,file='zgoubi_searchDA-In.dat')

      jo = 1
 2    continue

        write(icho,*) ' Now building zgoubi.dat for jo=',jo
        call flush2(icho,.false.)

        if( (xplus(jo)-xminu(jo)).le.prec ) then
          write(icho,fmt='(a,1x,i5,a)')  'Particle #',jo,
     >                                'rejected, has zero DA'
          call flush2(icho,.false.)
          write(lunSto,fmt='(a,1x,i5,a)')  '# Particle number',jo,
     >                                'rejected, has zero DA'
          if(jo .ge. nTraj) goto 60
          jo = jo+1
          goto 2 
        else
          write(icho,fmt='(a,1x,i5,2a)')  'Particle #',jo,' ,',
     >              'DA search started'
          call flush2(icho,.false.)
        endif

        open(unit=lunW,file='zgoubi.dat')
        rewind(lunW)
        rewind(lunR)

C Read till "KOBJ"
        do i=1,4
          read(lunR,fmt='(a)') txt132
          write(lunW,*) txt132
          write(*,*) txt132
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
C Complete OBJET with the line with '1'
 37     read(lunR,fmt='(a)',end=62) txt132
        txt132 = txt132(debstr(txt132):finstr(txt132))
        if(txt132(1:2) .eq. '1 ') then
          goto 37
        else
          backspace(lunR)
        endif
        write(lunW,*) ' 1 '
        txt132 = '''FAISTORE'''
        write(lunW,*) txt132
        txt132 = 'b_zgoubi.fai'
        write(lunW,*) txt132
        txt132 = ' 1'
        write(lunW,*) txt132
C Completes zgoubi.dat with the rest of zgoubi_searchDA-In.dat
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
            endif
          if    (strcon(txt132,'''REBELOTE''',10,
     >                                        IS) ) then 
            write(lunW,*) txt132   
            read(lunR,fmt='(a)',end=62) txt132
            write(txt132(1:6),fmt='(i5)') nTurn
            txt132 = txt132(1:6)//'  0.0  99'
            write(lunW,*) txt132
            txt132 = '''END'''
            write(lunW,*) txt132
            goto 10
          endif
          if    (strcon(txt132,'''END''',5,
     >                                     IS) ) then 
              txt132 = '''REBELOTE'''
              write(lunW,*) txt132
              write(txt6,fmt='(I6)') nTurn
              txt132 = txt6//'  0.0  99'
              write(lunW,*) txt132
              txt132 = '''END'''
              write(lunW,*) txt132
              goto 10
          endif

          write(lunW,*) txt132   
          write(*,*) txt132   

        goto 1    !!go-on completing zgoubi.dat with the rest of zgoubi_searchDA-In.dat

 10     continue
        close(lunW)

C Scan x ; 
C scan from x=0 to max>0 and then from x= -dx to min<0
        kdx = 0
        idx = 1
        klost = 0   ! used to allow zero DA over klost successive x values, before giving up 
        dx = (xplus(jo) - xminu(jo))/float(2*Ndx+1)
        isign = +1
        if( xplus(jo) .lt. (xst0+prec) ) then 
          xst0 = xplus(jo) - dx
          write(icho,*)  'xplus(jo).lt.xst0+prec)', 
     >              'hence  xst0 = xplus(jo) - dx'
          call flush2(icho,.false.)
        else
          write(*,*) ' Scan x **** ',jo, dx ,xplus(jo), xst0
          xst0 = x(jo)
          write(icho,*)  'xplus(jo).ge.xst0+prec)', 
     >              'hence  xst0 = x(jo)'
          call flush2(icho,.false.)
c          isign = -1
c          write(icho,fmt='(a,1x,i5)')  'Particle #',jo,
c     >    'launched with sole negative scan, no DA beyond x0. '
c          write(lunSto,fmt='(a,1x,i5)')  '# Particle number ',jo,
c     >    'launched with sole negative scan, zero positive region. '
        endif
 21     continue

        write(icho,fmt='(a,1x,i5,a,1x,1p,e14.6,1x,2i8,a
     >   ,/,6(1x,e12.4))') 
     >    ' Now looking for max Z of particle #',jo
     >   ,'  dx, kdx, isign : ',dx, kdx, isign 
     >   ,' x(jo),xp(jo),z(jo),zp(jo),s(jo),d(jo) :'
     >   ,  x(jo),xp(jo),z(jo),zp(jo),s(jo),d(jo)
        call flush2(icho,.false.)

        xst = xst0 + kdx * dx
        xpst = xp(jo)
            zst = z(jo)
            du = zDA
            if(zDA.lt.1.) du = duStrt(jo)
        if(kdx.eq.0
     >    .or. chs) then
            du = duStrt(jo)  ! cm
        endif
        zpst = zp(jo)

        write(icho,fmt='(a,1x,4(1x,e14.6))')
     >   ' first Z-stab limit search for xst, xpst, zst, zpst : ',
     >    xst, xpst, zst, zpst 
        call flush2(icho,.false.)

        call stabLim(HV,prec,jo,xst,xpst,zst,zpst,s(jo),d(jo),du,
     >                  zDA,zpDA,dDA,kex,ok)

C Here, xst0 is assumed to be the closed orbit coordinate
      storb(1,idx) = xst-xst0
      storb(2,idx) = zDA
      storb(3,idx) = dDA
      storb(4,idx) = jo
      storb(5,idx) = xst0      
      write(lunSto,*) '# ',xst-xst0,zDA,dDA,xst0,jo,klost
      call flush2(lunSto,.false.)

      chs = .false.
      if(zDA.le.prec) then
        klost = klost + 1
        if(klost.gt.klostM) then 
          klost = 0
          if(isign.eq.-1) then
            write(*,*) ' End run  jo / nTraj :    ', jo,' / ',nTraj
            write(icho,*) ' End run  jo / nTraj :    ', jo,' / ',nTraj
            call flush2(icho,.false.)
            kdxm = -kdx
C Ends run for particle jo. Prviously filled 'storb' will be printed 
C ordered from xmin to xmax. 
            call imp(storb,kdxp,kdxm,lunSto)
C            call imp(storb,kdx,lunSto)
            if(jo .ge. nTraj) goto 60
            jo = jo+1
            goto 2   
          else
            chs = .true.
            kdxp = kdx+1 
            isign = -1
            kdx = 0
          endif
        endif
      else
        klost = 0
      endif
      
      idx = idx+1
      kdx = kdx + isign
      goto 21

  62  continue
      write(6,*) ' Stopped sarchDA upon EOF in READ unit ...'
      write(6,*) ' Suggestion : check that zgoubi_searchDA-In.dat'
      write(6,*) '             is well completed...'
      goto 98
  60  continue
      close(lunR)
      close(lunW)

 999  continue
      write(6,*) ' *** Error upon reading lunData '
      goto 99

 998  continue
      write(6,*) ' *** End of file lunData reached'
      goto 99

 99   continue
      write(6,*) ' End of execution upon xst = final value'
 98   continue
      close(lunSto)
      close(icho)
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
      subroutine stabLim(HV,prec,jo,xst,xpst,zco,zpco,s,d,du,
     >                      zst,zpst,dst,kex,ok)
      character*(*) HV
      logical lost, ok, first
      parameter (lunR=13,lunW=14)
      character txt132*132, let*1, cmnd*200
      data first / .true. /
      data let / 'i' /

      iter = 0
      ok = .false.
      lost = .false.

C du stands for dx or dz (cm)
c      du = 5.   !!cm 
c          if(HV .eq. 'V') then
            zco = zco + du
            zst = zco 
c            xst = xco 
c          elseif(HV .eq. 'H' .or. HV .eq. 'Hz') then
c            xco = xco + du
c            xst = xco 
c            zst = zco
c          endif
c      xpst = xpco
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
c          write(*,fmt='(1p,a,4e14.6,e9.1,e12.4)') 
c     >    'Initial object :', xst,xpst,zco,zpco,sst,dst
c     >    'Initial object :', xco,xpco,zco,zpco,sst,dst
c          write(*,fmt='(1p,a,4e14.6,e9.1,e12.4,4a,i4)') 
c     >    'New object :', xst,xpst,zst,zpst,sst,dst,' ','''','A','''',jo
C Completes zgoubi.dat with the rest of searchStabLim.temp2
 19       continue
          read(lunR,fmt='(a)',end=10) txt132
          write(lunW,fmt='(a)') txt132   
        goto 19 
 10     continue
        close(lunR)
        close(lunW)
C-----------------------------------------------------------------------

C Run zgoubi with single traj
        cmnd = '~/zgoubi/source/zgoubi'
        write(*,*) cmnd
        call system(cmnd)
    
C Get initial Traj coord from zgoubi.res
        call getRes(
     >              lost)

C      write(*,*) ' x,xp,z,zp,s,d,lost : ',xst,xpst,zst,zpst,sst,dst,lost

        if(lost) then
          write(6,*) 
          write(6,*) '-------------------------------------------- '
          write(6,*) '--  Particle # ',jo,' lost, back to previous '
          write(6,*) '    zco+du -> zco,  du=',du,'   x, z, iter :',
     >                                         xst,zst,iter
          write(6,*)

            zco = zco - du
            zst = zco 

c                write(*,*) ' HV, xst, zst :', hv, xst, zst

          du = du/2. 
          if (du.le.prec) then
            ok = .true.
            lost = .true.
            call system('mv -f b_zgoubi-temp.fai b_zgoubi.fai')
            goto 99
          endif
          lost = .false.
          call system('mv -f b_zgoubi.fai b_zgoubi-temp.fai')
        endif      

            zco = zco + du
            zst = zco 

        iter = iter+1
c        if(iter.le.10) then
           goto 2
c        else
c           ok = .false.
c           goto 99
c        endif

 99   continue
      write(6,*) ' ///////////////////////////////////////////'
      write(6,*) ' End of run of particle  #jo =',jo,'  : '
      write(6,*) ' du, x, z, lost, iter : ',du,xst,zst,lost,iter
      write(6,*) ' ///////////////////////////////////////////'
      call system('rm searchStabLim.temp2')
      return
      end
      subroutine getRes(
     >                  lost)
C Get initial Traj coord, and average orbit, from zgoubi.res
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
c        write(*,*) '****  Found  part. lost, with x, z =',x, z
c        write(6,*)  ' ' 
        lost = .true.
      goto 99

 18   continue
      lost = .false.
c      write(*,*) ' sbr getRes * End of read, no loss *   '
      goto 99

 19   continue
      lost = .false.
c      write(*,*) ' sbr getRes * End of readPU upon read error *   '

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
      subroutine imp(storb,kdxp,kdxm,lunSto)
      parameter (ixMax=400)
      dimension storb(7,ixMax), temp(7,ixMax)
      nx = kdxp+kdxm
      do i = 1, nx
        do j=1,5
          temp(j,i) = storb(j,i)
        enddo
      enddo
      do i = 1, kdxm
        do j=1,5
          storb(j,i) = temp(j,nx-i+1)
        enddo
      enddo
      do i = 1, kdxp
        do j=1,5
          storb(j,kdxm+i) = temp(j,i)
        enddo
      enddo
      do i = 1, nx
        write(lunSto,*) (storb(j,i),j=1,5)
      enddo
      call flush2(lunSto,.false.)
      return
      end

