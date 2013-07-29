      parameter (MXX =100, MXY =100)
      dimension qx(MXX,MXY), qy(MXX,MXY), surf(MXX,MXY)
      dimension qxMad(MXX*MXY), qyMad(MXX*MXY)
      dimension strng(6,MXX*MXY)
      parameter (lunR=7,lunw=8,icho=10,isav=12)
      character*132 txt132
C prec on z-limit, number of turns, starting value for du
C                               cm              cm 
      data prec, nt, duStrt / 0.002,   999,      1. / 
C
      open(unit=icho,file='calcDA.echo')
      open(unit=isav,file='calcDA.out')
      write(isav,fmt='(a)') 
     > '#  qx  qy  qxMad  qyMad  DA  x+  x-  y-  y_0  y+  m,n,mn'

C Read strengths Ks(m,n) from prior MAD run 
      open(unit=lunR,file='scanTunes_Ks.out')
      read(lunR,fmt='(a132)',err=2,end=11) txt132
      write(*,*) txt132
      mn = 0
      nr = 0
 2    continue
        nr = nr+1
        read(lunR,*,err=2,end=11) (strng(iqua,mn),iqua=1,6), 
     >      qxMad(mn),qyMad(mn)
            write(88,*)     qxMad(mn),qyMad(mn)


        mn = mn+1
        goto 2
 11     continue
      nerr = nr-mn
      write(*,*) ' ok, read strengths, mn=',mn,', read errors =',nerr
      mMa = 41
      nMa = 41
      if(mMa .gt. MXX) stop 'increse size MXX'
      if(nMa .gt. MXY) stop 'increse size MXY'

      m = 1
      n = 1
      mn = 1
 1    continue

Create zgoubi_StabLim-In.dat with updated quad strengths
      call system('cp zgoubi_calcDA-In.dat zgoubi_StabLim-In.dat')
      open(unit=lunR,file='zgoubi_StabLim-In.dat')
      open(unit=lunw,file='zgoubi_StabLim-In.dat_temp')
      call newDat(lunR,lunw,mxx,mxy,strng,mn)
      close(lunw)
      close(lunR)
      call system('mv zgoubi_StabLim-In.dat_temp zgoubi_StabLim-In.dat')

! Stage 1 :   get paraxial tunes (Qm,Qn)
!   1.a create zgoubi.dat
!   1.b create traj in zgoubi.dat with paraxial x0 and paraxial y0
      open(unit=lunR,file='zgoubi_StabLim-In.dat')
      open(unit=lunw,file='zgoubi.dat')
      do i = 1, 5
        read(lunR,fmt='(a132)') txt132
        write(lunw,fmt='(a132)') txt132
      enddo
      read(lunR,fmt='(a132)') txt132
      write(lunw,*) '0.0001 .0 0.0001 .0 0. 1. ''o'''
      call gtEnd(lunR,lunw,icho)
      close(lunR)
      close(lunw)
!   1.c run zgoubi
      call system('~/zgoubi/source/zgoubi')
!   1.d compute tunes
      call system('~/zgoubi/struct/tools/tunesFromFai/tunesFromFai')  
!   1.e save tunes
      open(unit=lunw,file='tunesFromFai.out')
      read(lunw,fmt='(a132)') txt132
      read(lunw,*) x,xp,qx(m,n),qy(m,n)
      close(lunw)    
      write(icho,fmt='(2(1x,f10.6),2(1x,f8.4),3(1x,i4),2(1x,f10.6),a)') 
     >  qx(m,n),qy(m,n),qx(m,n)-qxMad(mn),qy(m,n)-qyMad(mn),m,n,mn
     > ,qxMad(mn),qyMad(mn)
     > ,' qx(m,n), qy(m,n), qx-qx_mad, qy-qy_mad'
       call flush2(icho)      

! Stage 2 :  search lower stable limit x-
!   2.a create searchStabLim.data for Hz lower limit x- (including subliminal z)
      open(unit=lunw,file='searchStabLim.data')
      write(lunw,fmt='(a,1p,e14.6,1x,i6,a)') 'Hz ',prec,nt,' 0 1 -1. 0'
      write(lunw,*) '#H/Hz/V  prec  #turn  trck #sampl  duStrt  addPrxl'
      close(lunw)
!   2.b search x-. This will also create zgoubi_StabLim-Out.dat_Hz
      call system(
     >'~/zgoubi/struct/tools/searchStabLim/searchStabLim')
!   2.c Read x-
      open(unit=lunw,file='zgoubi_StabLim-Out.dat_Hz')
      do i = 1, 5
        read(lunw,fmt='(a132)') txt132
      enddo
      read(lunw,*) xminu
      close(lunw)
      call system('mv 
     >    zgoubi_StabLim-Out.dat_Hz zgoubi_StabLim-Out.dat_xm')

! Stage 3 : search upper stable limit x+
!   3.a create searchStabLim.data for Hz upper limit x+ (including subliminal z)
      open(unit=lunw,file='searchStabLim.data')
      write(lunw,fmt='(a,1p,e14.6,1x,i6,a)') 'Hz ',prec,nt,' 0 1 1. 0'
      write(lunw,*) '#H/Hz/V  prec  #turn  trck #sampl  duStrt  addPrxl'
      close(lunw)
!   3.b search x+. This will also create zgoubi_StabLim-Out.dat_Hz
      call system(
     >'~/zgoubi/struct/tools/searchStabLim/searchStabLim')
!   3.c Read x+
      open(unit=lunw,file='zgoubi_StabLim-Out.dat_Hz')
      do i = 1, 5
        read(lunw,fmt='(a132)') txt132
      enddo
      read(lunw,*) xplus
      close(lunw)
      call system('mv 
     >    zgoubi_StabLim-Out.dat_Hz zgoubi_StabLim-Out.dat_xp')

! Stage 4 :  search  z-limit at x=0
!   4.a create searchStabLim.data for V  limit z0
      open(unit=lunw,file='searchStabLim.data')
      write(lunw,fmt='(a,1p,e14.6,1x,i6,a)') 'V ',prec,nt,' 0 1 1. 0'
      write(lunw,*) '#H/Hz/V  prec  #turn  trck #sampl  duStrt  addPrxl'
      close(lunw)
!   4.b search z0. This will also create zgoubi_StabLim-Out.dat_V
      call system(
     >'~/zgoubi/struct/tools/searchStabLim/searchStabLim')
!   4.c Read z0
      open(unit=lunw,file='zgoubi_StabLim-Out.dat_V')
      do i = 1, 5
        read(lunw,fmt='(a132)') txt132
      enddo
      read(lunw,*) bid, bid, z0
      close(lunw)
      call system('mv 
     >  zgoubi_StabLim-Out.dat_V zgoubi_StabLim-Out.dat_z0')

      write(icho,*) ' x+, x-, z(0) : ',xplus, xminu, z0, m, n

! Stage 5 :
!   5.a create zgoubi_StabLim-In.dat & searchStabLim.data for V  limit at x-/2, zm
      open(unit=lunR,file='zgoubi_StabLim-In.dat')
      open(unit=lunw,file='zgoubi.dat')
      do i = 1, 5
        read(lunR,fmt='(a132)') txt132
        write(lunw,fmt='(a132)') txt132
      enddo
      read(lunR,fmt='(a132)') txt132
      write(lunw,fmt='(1p,e12.4,a)') xminu/2.,' .0 .0 .0 0. 1. ''o'''
      call gtEnd(lunR,lunw,icho)
      close(lunR)
      close(lunw)
      open(unit=lunw,file='searchStabLim.data')
      write(lunw,fmt='(a,1p,e14.6,1x,i6,a)') 'V ',prec,nt,' 0 1 1. 0'
      write(lunw,*) '#H/Hz/V  prec  #turn  trck #sampl  duStrt  addPrxl'
      close(lunw)
      call system('mv zgoubi.dat zgoubi_StabLim-In.dat')
!   5.b search zm. This will also create zgoubi_StabLim-Out.dat_V
      call system(
     >'~/zgoubi/struct/tools/searchStabLim/searchStabLim')
!   5.c Read zm
      open(unit=lunw,file='zgoubi_StabLim-Out.dat_V')
      do i = 1, 5
        read(lunw,fmt='(a132)') txt132
      enddo
      read(lunw,*) bid, bid, zm
      close(lunw)
      call system('mv 
     >  zgoubi_StabLim-Out.dat_V zgoubi_StabLim-Out.dat_zm')

! Stage 6 : 
!   6.a create zgoubi_StabLim-In.dat & searchStabLim.data for V  limit at x+/2, zp
      open(unit=lunR,file='zgoubi_StabLim-In.dat')
      open(unit=lunw,file='zgoubi.dat')
      do i = 1, 5
        read(lunR,fmt='(a132)') txt132
        write(lunw,fmt='(a132)') txt132
      enddo
      read(lunR,fmt='(a132)') txt132
      write(lunw,fmt='(1p,e12.4,a)') xplus/2.,' .0 .0 .0 0. 1. ''o'''
      call gtEnd(lunR,lunw,icho)
      close(lunR)
      close(lunw)
      open(unit=lunw,file='searchStabLim.data')
      write(lunw,*) 'V  ',prec,nt,'  0  1  1.  0'
      write(lunw,*) '#H/Hz/V  prec  #turn  trck #sampl  duStrt  addPrxl'
      close(lunw)
      call system('mv zgoubi.dat zgoubi_StabLim-In.dat')
!   6.b search zp. This will also create zgoubi_StabLim-Out.dat_V
      call system(
     >'~/zgoubi/struct/tools/searchStabLim/searchStabLim')
!   6.c Read zp
      open(unit=lunw,file='zgoubi_StabLim-Out.dat_V')
      do i = 1, 5
        read(lunw,fmt='(a132)') txt132
      enddo
      read(lunw,*) bid, bid, zp
      close(lunw)
      call system('mv 
     >  zgoubi_StabLim-Out.dat_V zgoubi_StabLim-Out.dat_zp')

! Now compute DA for (Qm,Qn)
      surf(m,n) = xplus /2. * ((z0+zp)/2. + zp/2.) +  
     >          (-xminu)/2. * ((z0+zm)/2. + zm/2.) 

      write(isav,fmt='(1p,10e12.4,3(1x,i4))') 
     > qx(m,n), qy(m,n), qxMad(mn), qyMad(mn), surf(m,n), 
     > xplus,xminu,zm,z0,zp, m,n,mn
      call flush2(isav)

      mn = mn+1
      m = m+1
      if(m.gt.mMa) then
        if(n.eq.nMa) goto 10
        m = 1
        n = n+1
      endif

      goto 1

 10   continue
c      write(isav,*) '///////////'
c      write(isav,fmt='(1p,3e14.6,2(1x,i4))') 
c     >    ((qx(m,n), qy(m,n), surf(m,n), m,n,m=1,mMa),n=1,nMa)
      close(isav)
      write(*,*) '  '
      write(*,*) ' End of job... '
      write(*,*) '  '

C Restore original file
      call system('cp zgoubi_StabLim-In.dat_temp zgoubi_StabLim-In.dat')
      stop
      end
      subroutine gtEnd(lunR,lunw,icho)
      character*132 txt132
 1    continue
        read(lunR,fmt='(a132)',err=98,end=99) txt132
        write(lunw,fmt='(a132)') txt132
        goto 1
 99   continue
      write(*,*) ' End reached in zgoubi.dat'
      return
 98   continue
      write(*,*) ' Read error in zgoubi.dat'
      return
      end
      subroutine newDat(lunR,lunw,mxx,mxy,strng,mn)
      dimension strng(6,mxx*mxy)
      character*132 txt132
      character*20 tx1, tx2, tx3, tx4, tx5, tx6
      logical strcon
      parameter (tx1='QUAD      QFI')
      parameter (tx2='QUAD      QDI')
      parameter (tx3='QUAD      QDMS1')
      parameter (tx4='QUAD      QFMS1')
      parameter (tx5='QUAD      QDMS2')
      parameter (tx6='QUAD      QFMS2')
      parameter (zro=0.d0)
      integer debstr, finstr
      
C Read Brho (assuming OBJET/KOBJ=2/IMAX=1)
      read(lunR,err=999,end=998,fmt='(a)') txt132
      txt132 = txt132(debstr(txt132):finstr(txt132))
      write(lunw,*) txt132
      read(lunR,err=999,end=998,fmt='(a)') txt132
      txt132 = txt132(debstr(txt132):finstr(txt132))
      write(lunw,*) txt132
      read(lunR,*,err=999,end=998) Brho
      write(lunw,*) Brho
      read(lunR,err=999,end=998,fmt='(a)') txt132
      txt132 = txt132(debstr(txt132):finstr(txt132))
      write(lunw,*) txt132
      read(lunR,err=999,end=998,fmt='(a)') txt132
      txt132 = txt132(debstr(txt132):finstr(txt132))
      write(lunw,*) txt132
      read(lunR,err=999,end=998,fmt='(a)') txt132
      txt132 = txt132(debstr(txt132):finstr(txt132))
      write(lunw,*) txt132
      read(lunR,err=999,end=998,fmt='(a)') txt132
      txt132 = txt132(debstr(txt132):finstr(txt132))
      write(lunw,*) txt132
C End of OBJET

      ic1 = 0
      ic2 = 0
      ic3 = 0
      ic4 = 0
      ic5 = 0
      ic6 = 0

 1    continue
        read(lunR,err=999,end=998,fmt='(a)') txt132
        txt132 = txt132(debstr(txt132):finstr(txt132))
        write(lunw,*) txt132
        if    (STRCON(txt132,TX1,
     >                           is)) then
          read(lunR,err=999,end=998,fmt='(a132)') txt132
          txt132 = txt132(debstr(txt132):finstr(txt132))
          write(lunw,*) txt132
          read(lunR,*,err=999,end=998) xl, rad, b0, b1
          iq = 1
          b1 = strng(iq,mn)*Brho*rad/1e4  
          write(lunw,fmt='(2(f10.4,1x),f7.1,1x,f12.6,8(1x,f2.0),a4)') 
     >              xl,rad,b0,b1,zro,zro,zro,zro,zro,zro,zro,zro,' new'
          ic1 = ic1+1

        elseif(STRCON(txt132,TX2,
     >                           is)) then
          read(lunR,err=999,end=998,fmt='(a132)') txt132
          txt132 = txt132(debstr(txt132):finstr(txt132))
          write(lunw,*) txt132
          read(lunR,*,err=999,end=998) xl, rad, b0, b1
          iq = 2
          b1 = strng(iq,mn)*Brho*rad/1e4  
          write(lunw,fmt='(2(f10.4,1x),f7.1,1x,f12.6,8(1x,f2.0),a4)') 
     >              xl,rad,b0,b1,zro,zro,zro,zro,zro,zro,zro,zro,' new'
          ic2 = ic2+1

        elseif(STRCON(txt132,TX3,
     >                           is)) then
          read(lunR,err=999,end=998,fmt='(a132)') txt132
          txt132 = txt132(debstr(txt132):finstr(txt132))
          write(lunw,*) txt132
          read(lunR,*,err=999,end=998) xl, rad, b0, b1
          iq = 3
          b1 = strng(iq,mn)*Brho*rad/1e4  
          write(lunw,fmt='(2(f10.4,1x),f7.1,1x,f12.6,8(1x,f2.0),a4)') 
     >              xl,rad,b0,b1,zro,zro,zro,zro,zro,zro,zro,zro,' new'
          ic3 = ic3+1

        elseif(STRCON(txt132,TX4,
     >                           is)) then
          read(lunR,err=999,end=998,fmt='(a132)') txt132
          txt132 = txt132(debstr(txt132):finstr(txt132))
          write(lunw,*) txt132
          read(lunR,*,err=999,end=998) xl, rad, b0, b1
          iq = 4
          b1 = strng(iq,mn)*Brho*rad/1e4  
          write(lunw,fmt='(2(f10.4,1x),f7.1,1x,f12.6,8(1x,f2.0),a4)') 
     >              xl,rad,b0,b1,zro,zro,zro,zro,zro,zro,zro,zro,' new'
          ic4 = ic4+1

        elseif(STRCON(txt132,TX5,
     >                           is)) then
          read(lunR,err=999,end=998,fmt='(a132)') txt132
          txt132 = txt132(debstr(txt132):finstr(txt132))
          write(lunw,*) txt132
          read(lunR,*,err=999,end=998) xl, rad, b0, b1
          iq = 5
          b1 = strng(iq,mn)*Brho*rad/1e4  
          write(lunw,fmt='(2(f10.4,1x),f7.1,1x,f12.6,8(1x,f2.0),a4)') 
     >              xl,rad,b0,b1,zro,zro,zro,zro,zro,zro,zro,zro,' new'
          ic5 = ic5+1

        elseif(STRCON(txt132,TX6,
     >                           is)) then
          read(lunR,err=999,end=998,fmt='(a132)') txt132
          txt132 = txt132(debstr(txt132):finstr(txt132))
          write(lunw,*) txt132
          read(lunR,*,err=999,end=998) xl, rad, b0, b1
          iq = 6
          b1 = strng(iq,mn)*Brho*rad/1e4  
          write(lunw,fmt='(2(f10.4,1x),f7.1,1x,f12.6,8(1x,f2.0),a4)') 
     >              xl,rad,b0,b1,zro,zro,zro,zro,zro,zro,zro,zro,' new'
          ic6 = ic6+1

        endif

        goto 1

 997  continue
      write(*,*) ' newDat :  end of job... '
      goto 99

 999  continue
      write(*,*) ' *** newDat : error upon reading  '
      goto 99

 998  continue
      write(*,*) ' *** newDat : end of file reached'
      goto 99

 99   continue
      write(*,*) '# of quads changed, 1 to 9 : ',ic1,ic2,ic3,ic4,ic5,ic6
      return
      end
      FUNCTION STRCON(STR,STR2,
     >                         IS)
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
1     CONTINUE
        DEBSTR=DEBSTR+1
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
      SUBROUTINE FLUSH2(IUNIT)
      CHARACTER TXT80
      IF(IUNIT.EQ.6.OR.IUNIT.EQ.5) RETURN 
      BACKSPACE(IUNIT)
        READ(IUNIT,FMT='(A80)') TXT80
      RETURN
      END
