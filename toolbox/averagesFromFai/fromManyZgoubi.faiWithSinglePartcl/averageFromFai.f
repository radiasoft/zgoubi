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

      CHARACTER(30) namFai, namfic, namfpx, folder

C Number of variables (the cols in .fai) to be, each independently, averaged over 
C the mxfai .fai files
      parameter (ncol=10)
      dimension icol(ncol)
      dimension avrg(ncol,mxpss), sig(ncol,mxpss)
      character(200) txt200
      logical strcon, binary

C      data namFai / 'b_zgoubi.fai'  /
C      data binary / .true. /
      data namFai / 'zgoubi.fai'  /
      data binary / .false. /
       data ipmax / -9999 /
           save ipmax

      write(*,*) '--------------------------------------------'
      write(*,*) 'Compute average of each col(i) at each pass,'//
     >' over several .fai files, a single particle each.'

      call system('rm -f fort.*')
      call system('mv -f averageFromFai.out averageFromFai.out_old')

      if(idluni(lunW)) 
     >  open(unit=lunW,file='averageFromFai.Out')

      if(idluni(lunRR)) 
     >  open(unit=lunRR,file='averageFromFai.in')
C folder can be '.\', or 'Run*', etc.
C namfai = something.fai for instance
C col(i=1,ncol) can be 10, 11 for Y, T, 14 for s, 
C 20, 21, 22 for SX, SY, SZ, etc. 
      read(lunRR,*,err=99,end=99) 
     >   folder,namfai,nbcol,(icol(i),i=1,nbcol)
      close(lunRR)
      write(*,*) '# ', trim(folder),' ',trim(namfai),' ',
     >  nbcol,(icol(i),i=1,nbcol),
     >   '  folder,namfai,ncol,(icol(i),i=1,nbcol)'
      write(lunW,*) '# ', trim(folder),' ',trim(namfai),' ',
     >  nbcol,(icol(i),i=1,nbcol),
     >   '  folder,namfai,ncol,(icol(i),i=1,nbcol)'

      avrg = 0.d0
      sig = 0.d0

      call system('ls -l . | cat > listFiles')

      if(idluni(nlst)) 
     >open(unit=nlst,file='listFiles')

      nbf = 0
      iout = 0
 1    continue
C      namfai = 'zgoubi.fai'
      read(nlst,fmt='(a)',err=10,end=10) txt200
c        write(6,*) trim(txt200)
c        write(*,*) trim(folder)
c        write(*,*) ' strcon : ',strcon(txt200,trim(folder),
c     >                             is),is 
c           read(*,*)
      if(strcon(txt200,folder,
     >                        is)) then
c        write(*,*) ' trim(txt200) ',trim(txt200)
c        write(*,*) ' trim(txt200(is:200)) ',trim(txt200(is:200))
c        write(*,*) ' folder **',folder
c        write(*,*) '***',trim(txt200(is:200)),'***',trim(namfai)
        write(*,*) trim(txt200(is:200))//'/'//trim(namfai)
     >  ,' - file #',nbf
cc               read(*,*)
        if(idluni(lunR)) open(unit=lunR,file=
     >  trim(txt200(is:200))//'/'//trim(namfai),err=10)

c        write(*,*) ' Now opened ',
c     >  trim(txt200(is:200))//'/'//trim(namfai)

c        write(*,*) 'Now getting kex < 0 - if any. File is #',nbf
        CALL HEADER(lunR,6,4,binary,*98)
        call getkex(lunR,
     >                  iex,ip,ipmax)
        kex(nbf+1) = iex
        if(iex .lt. 0) then
              iout = iout + 1
        else
          rewind(lunR)
          CALL HEADER(lunR,6,4,binary,*98)
          call sumUp(lunR,nbcol,icol,kex,
     >                            avrg,sig,ipmax,npart)
        endif
        nbf = nbf + 1
        write(*,*) 'Numb of prtcls out / # of files : ',iout,'/',nbf
     >  ,' Pass# / max. pass#  ip / ipmax : ',ip,' / ',ipmax
     >  ,',     kex = ',iex
c           if(nbf.eq.10) goto 10
c        goto 1
      endif
        close(lunR)
      goto 1

 10   continue
      close(nlst)
      write(*,*) ' '
      write(*,*) ' nbcol = ',nbcol
      write(*,*) ' Done going through folder files'
      write(*,*) ' Grabbed ', nbf ,' .fai files'
      write(*,*) ' Nmber of praticles out : ',iout
      write(*,*) ' Max number of passes found, ipmax : ',ipmax
      write(*,*) ' '

C      write(lunW,*) ' ip, avrg_col1,sig_col1, avrg_col2,sig_col2, ...'
      do ip = 1, ipmax
        do ic = 1, nbcol
c           if(ic.eq.1) write(*,*) ' --- ',avrg(ic,ip),ip,icol(ic)
c            read(*,*)
           if(icol(ic) .eq. 38 .or. icol(ic) .eq. 25) then
           else
            avrg(ic,ip) = avrg(ic,ip)  /dble(nbf-iout)
            sig(ic,ip) = sqrt(sig(ic,ip)/dble(nbf-iout) -avrg(ic,ip)**2)
          endif
        enddo
        write(lunW,fmt='(i7,1x,1p,7(e12.4,1x),3(1x,i7))') 
     >  nint(avrg(1,ip)),avrg(2,ip)
     >  ,(avrg(ic,ip),sig(ic,ip),ic=3,nbcol),nbf,nbf-iout,ip
      enddo
c      nbcol = 5
c      do ip = 1, ipmax
c        write(*,*) ' lunW = ',lunW, ip
c           read(*,*)
c        write(lunW,fmt='(i5,1x,6(e12.4,1x),i7)') 
c     >  ip,(avrg(ic,ip),sig(ic,ip),ic=2,nbcol),nbf
c      enddo

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

      subroutine sumUp(lunR,nbcol,icol,iex,
     >                   avrg,sig,ipmax,npart)
      implicit double precision (a-h,o-z)
      parameter(mxpss=1000000)
      parameter (ncol=10)
      dimension avrg(ncol,mxpss), sig(ncol,mxpss),npart(*)
      dimension icol(*)
      parameter(mxfai=9999)
      dimension iex(mxfai)
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

      CHARACTER(30) namfic, namfpx

C Number of variables (the cols in .fai) to be, each independently, averaged over 
C the mxfai .fai files
      dimension coll(50)

C ip is ipass
        ipmax = -9999999
        ip = 1  
 11     continue

          READ(lunR,*,ERR=4,END=3) 
C     >    KEX, d,y,t
     >    KEX,(coll(i),i=2,38)

c          if(iex(ip) .lt. 0) goto 12

c     >    KEX,(FO(J),J=1,7),
c     >    (F(J),J=1,7), 
c     >    (SO(J),J=1,4),sx, sy, sz, sn,
c     >    ENEKI, ENERG, 
c     >    IT, IREP, SORT, AMQ1,AMQ2,AMQ3,AMQ4,AMQ5, RET, DPR, PS,
c     >    BORO, IPSS, NOEL ,KLEY,LBL1,LBL2,LET 

c           write(*,*) kex,coll
c             read(*,*)

          if(ip.gt.mxpss) then
            write(*,*) ' Job stopped upon exceeded turn #'
            goto 2
          endif

          if(kex .le. 0) then 
C Most probably particles will be rejected by the test below, "vjn/epsyrms.gt.cutOff"
C before reaching here...
            koff = koff + 1
            write(lunW,*) '# pass# ',ipss,',  ifai # ',ifai,',  ip = '
     >      ,ip,',  prtcl has kex<0 =>  Now open next .fai.'
            close(lunR)
            npart(ip) = npart(ip) - 1
            goto 3
          endif

          if(ip.gt.ipmax) then 
            ipmax=ip
            namfpx = namfic
          endif

          if(npart(ip).gt.nprtmx) then 
            nprtmx = npart(ip)
            ipssmx = ip
          endif

          do i = 1, nbcol
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


      subroutine getkex(lunR,
     >                       kex,ip,ipmax)
      implicit double precision (a-h,o-z)
      parameter(mxfai=9999)
      PARAMETER ( PI=4.d0*atan(1.d0) , DEUXPI=2.d0*PI )

C Number of variables (the cols in .fai) to be, each independently, averaged over 
C the mxfai .fai files
      dimension coll(50)

 11     continue
          READ(lunR,*,ERR=4,END=3) 
C     >    KEX
     >    KEX,(coll(i),i=2,37),ip
                  
            if(ip .gt. ipmax) ipmax = ip
             coll(38) = ip

               if(kex.lt.0) goto 3
        goto 11

 3      continue
          if(ip .lt. ipmax) kex = -99
      write(*,*) 'getkex stopped upon eof. kex : ',kex
         return

 4    continue
          if(ip .lt. ipmax) kex = -99
      write(*,*) 'getkex stopped upon read-error. kex : ',kex

         return
           end

