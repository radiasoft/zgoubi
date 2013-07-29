      implicit double precision (a-h,o-z)
      parameter(mxfai=9999)
      parameter(mxpss=1000000)
      dimension spin(4,mxpss), Gga(mxpss), npart(mxpss), ipass(mxpss)
      PARAMETER ( PI=3.1415926536 , DEUXPI=2.0*PI )
      logical idluni
      dimension sigEpsy(mxfai), vjmxn(mxfai)

      PARAMETER (mxj=7,MXS=4)
      DIMENSION FO(MXJ),F(MXJ),SO(MXS),SF(MXS)

      PARAMETER (LBLSIZ=8)
      CHARACTER*(LBLSIZ) LBL1, LBL2

      PARAMETER (KSIZ=10)
      CHARACTER KLEY*(KSIZ) 
      CHARACTER TX1*1
      CHARACTER*1 LET

      integer debstr, finstr
      character*5 txt5
      PARAMETER (lunR=7)

      CHARACTER*30 NAMFIC, namfpx
      data namfic / 'b_zgoubi.fai'  /

      data ifai1, ifai2 / 0, 4031  / 
C      data ifai1, ifai2 / 0, 48  / 

      write(*,*) '--------------------------------------------'
      write(*,*) 'Computation of <S_y> averaged over particles'

      call system('rm -f fort.*')
      call system('mv -f averageS.out averageS.out_old')

      if(idluni(lunW)) 
     >  open(unit=lunW,file='averageS.Out')

      if(idluni(lunRR)) 
     >  open(unit=lunRR,file='averageS.in')
      read(lunRR,*,err=99,end=99) 
     >        ifai1,ifai2,y0,yp0,alfy,bety,epsyrms,cutOff
      close(lunRR)
      write(*,*) '# '
     >   ,  ifai1,ifai2,y0,yp0,alfy,bety,epsyrms
     >  ,' ifai1, ifai2, y0, yp0, alfy, bety, epsyrms'
      write(*,*) '# '
     >    ,cutOff,'     cut-off/epsyrms'

      do ip = 1, mxpss
        spin(1,ip) = 0.d0
        spin(2,ip) = 0.d0
        spin(3,ip) = 0.d0
        spin(4,ip) = 0.d0
        npart(ip) = 0
      enddo

      do i = 1, mxfai
        vjmxn(i) = 0.d0
      enddo

      ifai = ifai1-1
      ipmax = -mxpss
      ioff = 0
      koff = 0
 1    continue

        ifai = ifai+1
        if(ifai.gt.ifai2) then
          write(*,*) ' Reached requested # of .fai files ;',
     >    ' now will compute average and save in averageS.out'
          goto 2
        endif
        if(ifai.gt.mxfai-1) then
          write(*,*) ' # of .fai files requested is too large ;',
     >    ' now will compute average and save in averageS.out'
          goto 2
        endif
        txt5 = ''
            if    (ifai.le.9) then
              write(txt5,fmt='(i1)') ifai
              txt5 = '0000'//txt5(DEBSTR(txt5):FINSTR(txt5))
            elseif(ifai.le.99) then
              write(txt5,fmt='(i2)') ifai
              txt5 = '000'//txt5(DEBSTR(txt5):FINSTR(txt5))
            elseif(ifai.le.999) then
              write(txt5,fmt='(i3)') ifai
              txt5 = '00'//txt5(DEBSTR(txt5):FINSTR(txt5))
            elseif(ifai.le.9999) then
              write(txt5,fmt='(i4)') ifai
              txt5 = '0'//txt5(DEBSTR(txt5):FINSTR(txt5))
            elseif(ifai.le.99999) then
              write(txt5,fmt='(i5)') ifai
              txt5 = txt5(DEBSTR(txt5):FINSTR(txt5))
            else
              stop ' Conflicting ifai <-> txt5 '
            endif
            namfic = 'b_zgoubi.fai'//'_'//txt5            
c            write(88,*) 'namfic, txt5, ifai: ', namfic,' ',txt5,' ',ifai
                
c        if(idluni(lunR)) then
          close(lunR)
          open(unit=lunR,file=namfic,
     >       FORM='UNFORMATTED',err=2)
          write(*,*) ' Now opened ',namfic
          if(ifai.eq.0) CALL HEADER(lunR,88,4,.TRUE.,*98)
c        else
c          stop ' PGM averageS : no idle unit'
c        endif
      
        ip = 0  
 11     continue
          READ(lunR,ERR=4,END=3) 
     >    KEX,(FO(J),J=1,7),
     >    (F(J),J=1,7), 
     >    (SO(J),J=1,4),sx, sy, sz, sn,
     >    ENEKI, ENERG, 
     >    IT, IREP, SORT, AMQ1,AMQ2,AMQ3,AMQ4,AMQ5, RET, DPR, PS,
     >    BORO, IPSS, NOEL ,KLEY,LBL1,LBL2,LET 

c          write(89,fmt='(1p,2(i5,1x),i7,4(1x,g13.4))') 
c     >    ifai,it,ip,
c     >    sx, sy, sz, sn
          ip = ip+1
          IPASS(ip) = ipss
          npart(ip) = npart(ip) + 1

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
C            npart(ip) = npart(ip) - 1
            goto 1
          endif

          if(ip.gt.ipmax) then 
            ipmax=ip
            namfpx = namfic
          endif

          if(npart(ip).gt.nprtmx) then 
            nprtmx = npart(ip)
            ipssmx = ip
          endif

          if(ipss.eq.1) then
            y = fo(4)  *1d-2 - y0
            yp = fo(5) *1d-3 - yp0
            vj0 = gmy*y*y + 2.d0*alfy*y*yp + bety*yp*yp
            sigEpsy(ifai) = vj0 
            ph = 0.d0
          endif

          y = f(4)  *1d-2 - y0
          yp = f(5) *1d-3 - yp0
          vjn = gmy*y*y + 2.d0*alfy*y*yp + bety*yp*yp
          sigEpsyn = vjn           
          if(vjmxn(ifai).lt.vjn) vjmxn(ifai)=vjn/ epsyrms

            Gga(ip) = ENERG /938.27203d0 * 1.7928474
            spin(1,ip) = spin(1,ip) +sx
            spin(2,ip) = spin(2,ip) +sy
            spin(3,ip) = spin(3,ip) +sz
            spin(4,ip) = spin(4,ip) +sn
      
C          if(vjn/epsyrms.gt.cutOff) then
          if(vj0/epsyrms.gt.cutOff) then
              ioff = ioff + 1
              write(lunW,*) '# pass# ',ipss,',  ifai # ',ifai,',  ip = '
     >        ,ip,',  prtcl reached cut-off. Now go to next one.'              
              write(89,fmt='(i9,6(1x,e15.6),3(1x,i4),2a)') ip 
     >        ,fo(4)*1d-2-y0, fo(5)*1d-3-yp0, vj0/epsyrms
     >        ,f(4)*1d-2-y0, f(5)*1d-3-yp0, vjn/epsyrms
     >        ,ioff,koff,ifai
     >        ,'  ip, yo-y0, ypo-yp0, vj0/epsyrms, yn-yo, ypn-ypo,'
     >        ,' vjn/epsyrms, ioff, koff, ifai'
              if(250*(ifai/250).eq.ifai) call flush2(89,.false.)
              goto 1
          endif

        goto 11

 3      continue
        write(*,*) ' End of current .fai reached, now go to next one'
        if(250*(ifai/250).eq.ifai .or. ifai.eq.ifai2) then 
          write(*,*) ' Intermediate averages saved in averageS.Out'
          write(lunW,*) ' '
          do ip = 1, ipmax
            write(lunW,fmt=
     >         '(1p,i7,2(1x,i6),5(1x,e17.8),2(1x,i5),1x,e13.4,2a)') 
     >      ip, ipass(ip), npart(ip)
     >      ,spin(1,ip)/dble(npart(ip)),spin(2,ip)/dble(npart(ip))
     >      ,spin(3,ip)/dble(npart(ip)),spin(4,ip)/dble(npart(ip))
     >      ,Gga(ip), ifai, ioff, cutoff
     >    ,' ip, ipass(ip), npart(ip), sX,sY,sZ,|s|, G.gma, ifai, ioff,'
     >      ,' cut-off' 
          enddo
        endif

      goto 1
 
 4    continue
      write(*,*) ' Job stopped upon eof or read-error'
 2    continue
      close(lunR)

c      do ip = 1, ipmax
c        spin(1,ip) = spin(1,ip)/dble(npart(ip))
c        spin(2,ip) = spin(2,ip)/dble(npart(ip))
c        spin(3,ip) = spin(3,ip)/dble(npart(ip))
c        spin(4,ip) = spin(4,ip)/dble(npart(ip))
c      enddo


      write(*,*) '------------------------------------------------'
      write(*,*) 'Total number of turns : ',ipass(ipmax)
      write(*,*) 'Max # of prtcls, ',nprtmx,
     >                     ', was reached at pass # ',ipssmx
      write(*,*) '# of prtcls with emittance > cutOff is ',ioff
      write(*,*) '# of prtcls with kex<0 is ',koff
      write(*,*) '------------------------------------------------'

        do ip = 1, ipmax
            write(lunW,fmt=
     >         '(1p,i7,2(1x,i6),5(1x,e17.8),2(1x,i5),1x,e13.4,2a)') 
     >      ip, ipass(ip), npart(ip)
     >      ,spin(1,ip)/dble(npart(ip)),spin(2,ip)/dble(npart(ip))
     >      ,spin(3,ip)/dble(npart(ip)),spin(4,ip)/dble(npart(ip))
     >      ,Gga(ip), ifai, ioff,  cutoff
     >    ,' ip, ipass(ip), npart(ip), sX,sY,sZ,|s|, G.gma, ifai, ioff,'
     >      ,' cut-off' 
        enddo

 22   continue
      write(*,*) ' Job completed !'
      close(lunW)
      close(lunR)
      stop
 
 98   continue
      stop 'Error during read header in b_zgoubi.fai'
 99   continue
      stop 'Error during read in averageS.in'
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

      WRITE(NW,FMT='(10X,''File header  ('',I1,
     >'' lines) : '')') N

      IF(.NOT.BINARY) THEN
        READ(NL,FMT='(A80)',ERR=99,END=99) TXT80
        WRITE(NW,FMT='(A)') TXT80
        READ(NL,FMT='(A80)',ERR=99,END=99) TXT80
        WRITE(NW,FMT='(A)') TXT80
      ELSE
        READ(NL,ERR=99,END=89) TXT80
        WRITE(NW,FMT='(A)') TXT80
        READ(NL,ERR=99,END=89) TXT80
        WRITE(NW,FMT='(A)') TXT80
      ENDIF
      IF(.NOT.BINARY) THEN
        DO 1 I=3, N
           READ(NL,FMT='(A)',ERR=99,END=89) TXT80
           WRITE(NW,FMT='(A)') TXT80
 1      CONTINUE
      ELSE
        DO 2 I=3, N
           READ(NL,ERR=99,END=89) TXT80
           WRITE(NW,FMT='(A)') TXT80
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
