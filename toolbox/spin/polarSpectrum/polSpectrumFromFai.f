Compute avergae <p> and <p^2>, to deduce average tune and tune spread from chroma
      implicit double precision (a-h,o-z)
      character(3) txti
      character(800) txt800

      parameter (mFldr=1024)
      parameter (c=2.99792458d8, am = 0.511d0, gyro = 1.159652E-03)
      dimension nbPass(mFldr)
      logical first
      integer debstr, finstr

      logical exs
      
      data first / .true. /
      data mxpass / 999999 /

      temp = mxpass
      inquire(file='polSpectrumFromFai.in',exist=exs)
      if(exs) open(unit=3, file='polSpectrumFromFai.in')
      read(3,*,err=22,end=22) mxpass
      goto 23
 22   continue
      mxpass = temp
 23   continue
      close(unit=3)      
      if(exs) then
        write(*,*) 'Pass # found in polSpectrumFromFai.in, ok !' 
      else
        write(*,*) 'Can specify pass # in polSpectrumFromFai.in...'
      endif
      write(*,*) 'Average is computed at pass # ',mxpass
      
      open(unit=2, file='polSpectrumFromFai.out')
      write(2,fmt='(a)') 
     >'# 1-folder#, 2-<a.gam>, 3-sig(a.gam), 4-pRef, 5-a.gam_ref, '
     >//'6-8-SX-Z, 9-|S|, 10-nPass(folder), 11-min(a.gma),'
     >//'  12-max(a.gma),  13-am,  14-a, nPass(folder)'

 88   continue
      write(*,*) 'Give number of folders ( >0 & <',mFldr,') :'
      read(*,*,err=88,end=88) nFldr
      if(nFldr .gt. mFldr) goto 88
      write(*,*) 'Give x and y chromas (normally 1., 1.) :'
      read(*,*,err=88,end=88) xix, xiy
 
      do i = 0, nFldr-1
 
C          write(txti,fmt='(i3.3)') i    
          write(txti,fmt='(i0)') i

          open(unit=1,file='Run'//trim(txti)//'/zgoubi.fai',
     >    action='read',err=77)
          write(*,*) 'Opened '//'Run'//trim(txti)//'/zgoubi.fai.'
     >    //' Now computing <p> and <p^2>. '
          nbPass(i+1) = 0
          maxPass = -9999
          Dav = 0.d0
          D2av = 0.d0
          Dmi = 1d10
          Dma = -1d10
       
C Read header
          do j = 1, 4
            read(1,fmt='(a)',err=10,end=10) txt800
          enddo          

 1        continue
            read(1,fmt='(a)',err=10,end=10) txt800

            read(txt800(125:151),*,err=10,end=10) dpp     ! = (p-p0)/p0
            pp = 1.d0 + dpp
            Dav = Dav + pp
            D2av = D2av + pp*pp
            if(pp .lt. Dmi) Dmi=pp
            if(pp .gt. Dma) Dma=pp

            read(txt800(647:655),*,err=10,end=10) ipass

            if(first) read(txt800(629:646),*,err=10,end=10) BORO
              first=.false.

            read(txt800(364:428),*,err=10,end=10) sx, sy, sz, sm

            SXf=sx; SYf=sy; SZf=sz; Smod=sm
            nbPass(i+1) = nbPass(i+1) +1
            if (ipass.gt. maxPass) maxPass = ipass

          goto 1

 10       continue          

          close(1)
          first = .true.
          pRef = BORO*(c/1d9)       ! MeV
          aGamRef =  gyro * pRef/am 
          avGa = Dav/dble(nbPass(i+1)) * aGamRef
          av2Ga = D2av/dble(nbPass(i+1)) * aGamRef* aGamRef
          sigGa = sqrt(av2Ga - avGa*avGa)
          Gami = Dmi * aGamRef
          Gama = Dma * aGamRef
          
          write(2,fmt='(i6,1x,1p,8(e14.6,1x),i8,1x,5(e14.6,1x))') 
     >    i+1, avGa, sigGa, pRef, aGamRef, SXf, SYf, SZf, Smod, 
     >    nbPass(i+1), Gami, Gama, BORO, am, gyro
          call flush(2)

 77       continue

      enddo



        stop
          end
      
      FUNCTION DEBSTR(STRING)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
