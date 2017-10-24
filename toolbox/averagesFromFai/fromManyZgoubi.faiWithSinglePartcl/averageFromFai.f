      implicit double precision (a-h,o-z)
      parameter(mxfai=9999)
      parameter(mxpss=1000000)
      dimension npart(mxpss), ipass(mxpss)
      PARAMETER ( PI=4.d0*atan(1.d0) , DEUXPI=2.d0*PI )
      logical idluni

      PARAMETER (mxj=7,MXS=4)
      DIMENSION FO(MXJ),F(MXJ), kex(mxfai)

      PARAMETER (LBLSIZ=8)
      CHARACTER*(LBLSIZ) LBL1, LBL2

      PARAMETER (KSIZ=10)
      CHARACTER KLEY*(KSIZ) 
      CHARACTER TX1*1
      CHARACTER*1 LET

      integer debstr, finstr

      CHARACTER(30) namFai, namfic, folder

C Number of variables (the cols in .fai) to be, each independently, averaged over 
C the mxfai .fai files
      parameter (ncol=10)
      dimension icol(ncol)
      dimension avrg(ncol,mxpss), sig(ncol,mxpss)
      character(200) txt200
      logical strcon, binary
      character(300) txfrm
      character(3) txncol
      character(20) namfol
      
C      data namFai / 'b_zgoubi.fai'  /
C      data binary / .true. /
      data namFai / 'zgoubi.fai'  /
      data binary / .false. /
       data ipmax / -9999 /
           save ipmax

      write(*,*) '--------------------------------------------'
      write(*,*) 'Compute average and sigma of each col(i) at each '//
     >'pass, over several .fai files, a single particle each '//
     >'(except for col=25=ENERG and col=38=PASS#, restituted as read).'

      call system('rm -f fort.*')
      call system('mv -f averageFromFai.out averageFromFai.out_old')

      if(idluni(lunW)) open(unit=lunW,file='averageFromFai.Out')
      if(idluni(lunW2)) open(unit=lunW2,file='averageFromFai.Out2')
      if(idluni(lunRR)) open(unit=lunRR,file='averageFromFai.in')
C folder can be '.\', or 'Run*', etc.
C namfai = something.fai for instance
C col(i=1,ncol) can be 10, 11 for Y, T, 14 for s, 
C 20, 21, 22 for SX, SY, SZ, etc. 
      read(lunRR,*,err=99,end=99) 
     >   folder,namfai,nbcol,(icol(i),i=1,nbcol),iprdmx
      close(lunRR)

c      write(*,*) '  iprdmx =',iprdmx
c       read(*,*)
      
      icol(nbcol+1) = 38   ! IPASS
      icol(nbcol+2) = 25   ! ENERG
      write(*,*) '#  Folder: ', trim(folder),'*  File: ',trim(namfai),
     >nbcol,(icol(i),i=1,nbcol),'   '//
     >'folder,namfai,ncol,(icol(i),i=1,nbcol) and 38 (IPASS), 25(ENERG)'
      write(lunW,*) '#  Folder: ', trim(folder),'*  File: ',trim(namfai)
     >,nbcol,(icol(i),i=1,nbcol),'   '//
     >'folder,namfai,ncol,(icol(i),i=1,nbcol) and 38 (IPASS), 25(ENERG)'
      nbcol2 = nbcol+2
      
      write(txncol,fmt='(i3.3)') 2*nbcol
      txfrm = '(1p,'//trim(txncol)//'(e12.4,1x)' 
      txfrm = trim(txfrm)//',i7,1x,e12.4'
      txfrm = trim(txfrm)//',3(1x,i7))'
      write(*,*) '**'//trim(txfrm)//'**'
              
      avrg = 0.d0
      sig = 0.d0

      call system('ls -l . | cat > listFiles')

      if(idluni(nlst)) open(unit=nlst,file='listFiles')

      nbf = 0
      iout = 0
 1    continue
        read(nlst,fmt='(a)',err=10,end=10) txt200
        if(strcon(txt200,folder,          ! work on Run# found next in list
     >                        is)) then
        read(txt200(is:200),*) namfol
        
        if(idluni(lunR)) open(unit=lunR,file=
     >    trim(txt200(is:200))//'/'//trim(namfai),err=10)

          CALL HEADER(lunR,6,4,binary,*98)
          call getiex(lunR,nbf,iprdmx,               ! check if particle out at some point in this Run#
     >                  iex,ip,ipmax)
          kex(nbf+1) = iex

c     write(*,*) ' iex = ',iex
c          read(*,*)

          if(iex .ne. 1) then
            iout = iout + 1
            write(*,*) namfol(1:8),
     >      '  iex .ne. 1 -> this particle(file# ',nbf,') not counted.'
     >      //' A total of iout/nbrOfFiiles = ',iout,'/',nbf+1
c            read(*,*)
         else
           rewind(lunR)
           CALL HEADER(lunR,6,4,binary,*98)
           call sumUp(lunR,nbcol2,icol,kex,
     >                           avrg,sig,ipmax,npart)
         endif

         if(100*(nbf/100) .eq. nbf)
     >        write(*,*) namfol(1:8),
     >      'Numb of prtcls out/# of files : ',iout,'/',nbf+1
     >   ,',    max. pass#/nb. stored : ',ip,' / ',ipmax
     >   ,',     iex = ',iex
         write(lunW2,*)namfol(1:8),
     >    'Numb of prtcls out/# of files : ',iout,'/',nbf+1
     >   ,',    max. pass#/nb. stored : ',ip,' / ',ipmax
     >   ,',     iex = ',iex
c          if(nbf.eq.10) goto 10
c        goto 1

         nbf = nbf + 1
         
      endif

      close(lunR)

      goto 1

 10   continue
      close(nlst)
      write(*,*) ' '
      write(*,*) ' nbcol = ',nbcol2-2
      write(*,*) ' Done going through folder files'
      write(*,*) ' Grabbed ', nbf ,' .fai files'
      write(*,*) ' Nmber of praticles out : ',iout
      write(*,*) ' Max number of passes found, ipmax : ',ipmax
      write(*,*) ' '

C     do ip = 1, ipmax
      ip = 1
      do while (ip .le. nint(avrg(nbcol+1,ip))) 
        do ic = 1, nbcol2
c           if(ic.eq.1) write(*,*) ' --- ',avrg(ic,ip),ip,icol(ic)
c            read(*,*)
           if(icol(ic) .eq. 38 .or. icol(ic) .eq. 25) then
           else
            avrg(ic,ip) = avrg(ic,ip)  /dble(nbf-iout)                     ! This is 
            sig(ic,ip) = sqrt(sig(ic,ip)/dble(nbf-iout) -avrg(ic,ip)**2)   ! in zgoubi units (cm, mrad, etc.)
          endif
        enddo

        write(lunW,fmt=txfrm) 
     >  (avrg(ic,ip),sig(ic,ip),ic=1,nbcol),
     >  nint(avrg(nbcol+1,ip)),avrg(nbcol+2,ip),     ! IPASS, ENERG
     >  nbf,nbf-iout,ip
        ip = ip + 1
      enddo

      write(*,*) ' Job completed !'
      write(*,*) ' '
      close(lunW)
      close(lunR)
      stop
 
 98   continue
      write(*,*) 'Error during read header in ',trim(namfic)
      stop 
 99   continue
      stop 'Error during read in averageFromFai.in'
      END
      FUNCTION IDLUNI(LN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      SUBROUTINE HEADER(NL,NW,N,BINARY,
     >                                 *)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL BINARY
      CHARACTER*80 TXT80

c      WRITE(NW,FMT='(10X,''File header  ('',I1,
c     >'' lines) : '')') N

      IF(.NOT.BINARY) THEN
        READ(NL,FMT='(A80)',ERR=99,END=99) TXT80
c        WRITE(NW,FMT='(A)') TXT80
        READ(NL,FMT='(A80)',ERR=99,END=99) TXT80
c        WRITE(NW,FMT='(A)') TXT80
      ELSE
        READ(NL,ERR=99,END=89) TXT80
c        WRITE(NW,FMT='(A)') TXT80
        READ(NL,ERR=99,END=89) TXT80
c        WRITE(NW,FMT='(A)') TXT80
      ENDIF
      IF(.NOT.BINARY) THEN
        DO 1 I=3, N
           READ(NL,FMT='(A)',ERR=99,END=89) TXT80
c           WRITE(NW,FMT='(A)') TXT80
 1      CONTINUE
      ELSE
        DO 2 I=3, N
           READ(NL,ERR=99,END=89) TXT80
c           WRITE(NW,FMT='(A)') TXT80
 2      CONTINUE
      ENDIF
      RETURN

 89   CONTINUE
      WRITE(6 ,*) 
     >'SBR HEADER : END of file reached while reading data file header'
      WRITE(NW,*) 
     >'SBR HEADER : END of file reached while reading data file header'
      RETURN 1

 99   CONTINUE
      WRITE(6,*) 
     >'*** SBR HEADER : READ-error while reading data file header'
      WRITE(6,*) '        ... Empty file ?'
      WRITE(NW,*) 
     >'*** SBR HEADER : READ-error while reading data file header'
      WRITE(NW,*) '        ... Empty file ?'
      RETURN 1
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
      SUBROUTINE FLUSH2(IUNIT,BINARY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      LOGICAL BINARY
      CHARACTER TXT80
      IF(IUNIT.EQ.6.OR.IUNIT.EQ.5) RETURN 
      BACKSPACE(IUNIT)
      IF(BINARY) THEN
        READ(IUNIT) TXT80
      ELSE
        READ(IUNIT,FMT='(A80)') TXT80
      ENDIF
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

      subroutine sumUp(lunR,nbcol2,icol,kex,
     >                   avrg,sig,ipmax,npart)
      implicit double precision (a-h,o-z)
      parameter(mxpss=1000000)
      parameter (ncol=10)
      dimension avrg(ncol,mxpss), sig(ncol,mxpss),npart(*)
      dimension icol(*)
      parameter(mxfai=9999)
      dimension kex(mxfai)
      dimension ipass(mxpss)
      PARAMETER ( PI=4.d0*atan(1.d0) , DEUXPI=2.d0*PI )
      logical idluni

      PARAMETER (mxj=7,MXS=4)
      DIMENSION FO(MXJ),F(MXJ),SO(MXS),SF(MXS)

      PARAMETER (LBLSIZ=8)
      CHARACTER*(LBLSIZ) LBL1, LBL2

      PARAMETER (KSIZ=10)
      CHARACTER KLEY*(KSIZ) 
      CHARACTER TX1*1
      CHARACTER*1 LET

      integer debstr, finstr

      CHARACTER(30) namfic

C Number of variables (the cols in .fai) to be, each independently, averaged over 
C the mxfai .fai files
      dimension coll(50)

C ip is ipass
C        ipmax = -9999999
        ip = 1  
 11     continue

          READ(lunR,*,ERR=4,END=3) iEX,(coll(i),i=2,38)

          if(ip.gt.mxpss) then
            write(*,*) ' Job stopped upon exceeded turn #'
            goto 2
          endif

          if(iex .ne. 1) then 
C Most probably particles will be rejected by the test below, "vjn/epsyrms.gt.cutOff"
C before reaching here...
               stop ' Pgm sumUp. Problem. This should never happen'
          endif

c          if(ip.gt.ipmax) then 
c            ipmax=ip
c            namfpx = namfic
c          endif

          if(npart(ip).gt.nprtmx) then 
            nprtmx = npart(ip)
            ipssmx = ip
          endif

          do i = 1, nbcol2
           if(icol(i) .eq. 38 .or. icol(i) .eq. 25) then
              avrg(i,ip) = coll(icol(i))
           else
              avrg(i,ip) = avrg(i,ip) + coll(icol(i))
              sig(i,ip) = sig(i,ip) + coll(icol(i))*coll(icol(i))
            endif
c              write(*,*) ' i = ',i,icol(i), coll(icol(i))
c              write(*,*) i, ip, avrg(i,ip), sig(i,ip)
c                     read(*,*)
          enddo

 12       continue
          ip = ip + 1
        goto 11

 3      continue
c        write(*,*) 
c     >' End of current .fai reached, now move on to next one'
         return

 4    continue
      write(*,*) ' Job stopped upon eof or read-error'
 2    continue

         return
           end

      subroutine getiex(lunR,nbf,iprdmx,
     >                       iex,ip,ipmax)
      implicit double precision (a-h,o-z)
      parameter(mxfai=9999)
      PARAMETER ( PI=4.d0*atan(1.d0) , DEUXPI=2.d0*PI )
      
C Number of variables (the cols in .fai) to be, each independently, averaged over 
C the mxfai .fai files
      dimension coll(50)

      ierr = -1
      
 11   continue
          READ(lunR,*,ERR=4,END=3) 
     >    IEX,(coll(i),i=2,37),ip
                  
            if(ip .gt. ipmax) ipmax = ip
             coll(38) = ip

             if(iex.ne.1) goto 3
             if(ip.ge.iprdmx) then
                ierr = 0
                goto 77
             endif
        goto 11

 3      continue
          ierr = 0
          if(ip .lt. ipmax) iex = -99
c      write(*,*) 'getiex stopped upon eof. iex : ',iex
         goto 77

 4    continue
        ierr = 4
      if(ip .lt. ipmax) iex = -99
c      write(*,*) 'getiex stopped upon read-error. iex : ',iex

 77   continue
c      write(*,*) ' getiex  ierr, ip, ipmax, iex, file# : ',
c     > ierr, ip, ipmax, iex, nbf+1 

c      write(*,*) ' getiex iprdmx ',iprdmx,iex
c       read(*,*)
      
 
      return
           end

