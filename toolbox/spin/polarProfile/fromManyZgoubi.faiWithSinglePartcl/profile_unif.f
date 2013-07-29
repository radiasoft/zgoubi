      implicit double precision (a-h,o-z)
      parameter(mxfai=9999)
      parameter(mxpss=1000000, mbin=100, mxph=40)
C      dimension spin(4,mbin,mxph,mxpss), npart(mxpss),ipass(mxpss)
      dimension spin(4,mbin), npart(mxpss),ipass(mxpss),nprt(mbin)
     > ,dx(mbin),xx(mbin)
      dimension iqij(mbin,mxfai)
      PARAMETER ( PI=3.1415926536 , DEUXPI=2.0*PI, Gp=1.79284735 )
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

      parameter (lunR=7)
      logical setvj0, first
      CHARACTER*30 NAMFIC, namfpx
      data namfic / 'b_zgoubi.fai'  /

C ifai1/2 : min/max xtensions of b_zgoubi.fai_...
      data ifai1, ifai2, nbin / 0, 4031, 40 / 
      data goalGg / 1.d10 / 
C      data ifai1, ifai2 / 0, 48  / 
      data setvj0, first / .false., .true. /

      call system('rm fort.*')
      call system('mv -f profileS.out profileS.out_old')

      write(*,*) '--------------------------------------------'
      write(*,*) 'Computation of poalrization y-profile'

      if(idluni(lunRR)) 
     >  open(unit=lunRR,file='profileS.in')
      read(lunRR,*,err=99,end=99) 
     >        ifai1,ifai2,y0,yp0,alfy,bety,epsyrms,cutoff,nbin
     >    ,goalGga,goalGgb,dGg
      goalGg = goalGga
      close(lunRR)
      if(nbin.gt.mbin) stop ' nbin is too large '

      N0 = ifai2 - ifai1
      gmy = (1.d0 + alfy*alfy)/bety      

      write(*,*) '# '
     >   ,  ifai1,ifai2,y0,yp0,alfy,bety,epsyrms,cutoff
     >  ,' ifai1, ifai2, y0, yp0, alfy, bety, epsyrms, cutoff'
      write(*,*) '# '
     >    ,nbin,goalGga,goalGgb,dGg
     >    ,' nbin, goalGga, goalGgb, dGg'

        frac = .2d0
        nbin =  5.d0 /frac + 1 

        dx(1) = frac*epsyrms  ! 4sigme of epsy is 98.17%
        xx(1) = dx(1)/2.d0
        dN = N0/dble(nbin)
        sumN = dN
C       write(89,*)i,xx(1)/epsyrms,dN, sumN,' i,xx(i)/epsyrms,dN, sumN'
      do i = 2, nbin
        dx(i) = dx(1)
        xx(i) = xx(i-1) + dx(i-1)/2.d0 + dx(i)/2.d0
        sumN = sumN + dN
C       write(89,*)i,xx(i)/epsyrms,dN, sumN,' i,xx(i)/epsyrms,dN, sumN'
      enddo

      if(idluni(lunW)) 
     >  open(unit=lunW,file='profileS.out')

C        do iph = 1, mxph
        do ij = 1, nbin
          spin(1,ij) = 0.d0
          spin(2,ij) = 0.d0
          spin(3,ij) = 0.d0
          spin(4,ij) = 0.d0 
          nprt(ij) = 0.d0 
          do ii = 1, mxfai
            iqij(ij,ii) = 0
          enddo
        enddo
C        enddo


      do ip = 1, mxpss
        npart(ip) = 0
      enddo

      first = .true.
      ifai = ifai1
      ipmax = -1
      ijmx = -1
      iout = 0
      koff = 0
      SZAvrg = 0.d0
      nbPrtcl = 0
 1    continue

        if(ifai.gt.ifai2) then
          write(*,*) ' Reached requested # of .fai files ;',
     >    ' now will save in profileS.out'
          goto 2
        endif
        if(ifai.gt.mxfai-1) then
          write(*,*) ' # of .fai files requested is too large ;',
     >    ' now will compute average and save in profileS.out'
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

          close(lunR)
          open(unit=lunR,file=namfic,
     >       FORM='UNFORMATTED',err=2)
          if(first) then
            CALL HEADER(lunR,88,4,.TRUE.,*98)
            first = .false.
          endif
      
        ip = 0  
 11     continue
          READ(lunR,ERR=4,END=31) 
     >    KEX,(FO(J),J=1,7),
     >    (F(J),J=1,7), 
     >    (SO(J),J=1,4),sx, sy, sz, sn,
     >    ENEKI, ENERG, 
     >    IT, IREP, SORT, AMQ1,AMQ2,AMQ3,AMQ4,AMQ5, RET, DPR, PS,
     >    BORO, IPSS, NOEL ,KLEY,LBL1,LBL2,LET 

          Gg = Gp * ENERG / AMQ1 
          if(Gg.gt.goalGg) goto 3

          if(ip.gt.mxpss) then
            write(*,*) ' Job stopped upon exceeded turn #'
            goto 2
          endif

          ip = ip+1
          IPASS(ip) = ipss
          npart(ip) = npart(ip) + 1

          if(kex .le. 0) then 
C Most probably particles will be rejected by the test below, "vjn/epsyrms.gt.cutOff"
C before reaching here...
            write(lunW,*) '# pass# ',ipss,',  ifai # ',ifai,',  ip = '
     >      ,ip,',  prtcl has kex<0 =>  Now open next .fai.'
            koff=koff+1
            write(*,*) 
     >      'Traj # ',it,',  kex<0  =>  now close '
     >      ,namfic(DEBSTR(namfic):FINSTR(namfic)),'. Total out : ',koff
            close(lunR)
            npart(ip) = npart(ip) - 1
            ifai = ifai+1
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
            setvj0 = .true.
c            write(89,fmt='(i4,5(1x,e15.6),1x,i4,a)') 
c     >       ip,fo(4)*1d-2,y0,fo(5)*1d-3,yp0,vj0,ifai,
c     >         '  ip,y,yp,vj0,ifai'
          endif

          y = f(4)  *1d-2 - y0
          yp = f(5) *1d-3 - yp0
          vjn = gmy*y*y + 2.d0*alfy*y*yp + bety*yp*yp
          sigEpsyn = vjn 
          
          if(vjmxn(ifai).lt.vjn) vjmxn(ifai)=vjn/ epsyrms

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
              goto 1
          endif

        goto 11

 3      continue

c        write(*,*) ' End of '
c     >  ,namfic(DEBSTR(namfic):FINSTR(namfic))        
c     >  ,' reached, now get next one'

c        write(*,*) ip,y,yp,vj0,ij,ijmx,'  ip,y,yp,vj0,ij,ijmx'

          do i = 1, nbin
            if(vj0.le.xx(i)) then
              ij = i
C          write(89,*) ifai,i,vj0,xx(i),' ifai,ij,vj0,xx(i)' 
              goto 33
            endif
          enddo
          ij = nbin

 33       continue

          if(ij.le.nbin) then
            if(ij.gt.ijmx) then 
              ijmx = ij
            endif
            iph = 1          
            spin(1,ij) = spin(1,ij) +sx
            spin(2,ij) = spin(2,ij) +sy
            spin(3,ij) = spin(3,ij) +sz
            spin(4,ij) = spin(4,ij) +sn
            nprt(ij) = nprt(ij) + 1
            iqij(ij,nprt(ij)) = ifai
            SZAvrg = SZAvrg  + sz
            nbPrtcl = nbPrtcl + 1            
          else
C             stop ' nbin exceed !'
            iout = iout+1
c            write(88,*) 'pass#, partcile#, ij = ',ip,ifai,ij
            write( *,*) 'pass#, partcile#, ij = ',ip,ifai,ij
          endif

          goto 32

 31   continue
      write(*,*) 'Exited READ upon EOF '

 32   continue
      if(Gg.lt.GoalGg) 
     >  write(*,*) 'WARNING : Ggamma, GoalGg = ',Gg, GoalGg


      ifai = ifai+1
      goto 1
 
 4    continue
      write(*,*) ' Job stopped upon eof or read-error'

 2    continue

      write(*,*) '------------------------------------------------'
      write(*,*) 'Total number of turns : ',ipmax
      write(*,*) '# of prtcls analyzed : ',nprtmx
      write(*,*) '# of prtcls with emittance > cutoff   is ',ioff
      write(*,*) '# of prtcls with kex<0 : ',koff
      write(*,*) '<SZ> = ',SZAvrg/dble(nbPrtcl),
     >                         '  with ',nbPrtcl,' partcls counted'
      write(*,*) '------------------------------------------------'

      max = ijmx
      Ntot = 0
      SZAvrg2 = 0.d0
      if(max.gt.nbin) max = nbin

        if(setvj0) then  
          do ij = 1, nbin
c             write(89,*) ij,nprt(ij),'  ij,nprt(ij)'
            if(nprt(ij) .ne.0) then
              Ntot = Ntot + nprt(ij)
              SZAvrg2 = SZAvrg2  + spin(3,ij)
              write(lunW,fmt='
     >        (1p,i3,6(1x,e17.8),3(1x,i6),1x,e17.8,1x,i8,1x,e17.8,2a)')
     >         ij, xx(ij)/epsyrms
     >        ,spin(1,ij)/dble(nprt(ij))
     >        ,spin(2,ij)/dble(nprt(ij))
     >        ,spin(3,ij)/dble(nprt(ij))
     >        ,spin(4,ij)/dble(nprt(ij)), SZAvrg2/dble(Ntot)
     >        ,ifai,  ip, nprt(ij),dx(ij),Ntot,Gg
     >        ,' ij, epsy/rms, sX, sY, sZ, |s|, SZAvrg2, ifai, ip,'
     >        ,' nprt(ij), dx, sum[dN], Ggamma' 
            endif
          enddo
        else
          stop ' Emittance unknown ! Job not completed'
        endif

      write(*,*) '------------------------------------------------'
      write(*,*) 'Second method : '
      write(*,*) '<SZ> = ',SZAvrg2/dble(Ntot),
     >                         '  with ',Ntot,' partcls counted'
      write(*,*) '------------------------------------------------'

      do ij = 1, nbin
        write(77,fmt='(a,i4,a)') 
     >     '# bin number ',ij,' contains following .fai :'
        do ii = 1, mxfai
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
          if(iqij(ij,ii).ne.0) 
     >       write(77,fmt='(3(i5,1x),1p,e13.4,2a)') ij,ii,
     >       xx(ij)/epsyrms, iqij(ij,ii),
     >       namfic(DEBSTR(namfic):FINSTR(namfic)),
     >       ' bin#, occurence#, epsy/sig_epsy, ifai value'
        enddo
      enddo
      close(77)
      call system('mv -f fort.77 profileS.out_iqij')

      write(*,*) ' Job completed !'

      close(lunW)
      close(lunR)

      stop
 
 98   continue
      stop 'Error during read header in b_zgoubi.fai'
 99   continue
      stop 'Error during read in profileS.in'
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
