c#ifdef __INTEL_COMPILER
c         use ifport
c#endif
      implicit double precision (a-h,o-z)
      character(800) txt800
      character(80) fileFai
      parameter (c=2.99792458d8, smu2s=1d-6)
      
      open(unit=1,file='runningAverage_fromFai.in')
      read(1,*,err=8,end=8) fileFai
      read(1,*,err=8,end=8) j1Pass
      goto 7
      close(1)
 8    continue
      fileFai = 'zgoubi.fai'
      j1Pass = 2000
 7    continue

      write(*,*) 'Read from file [pathTo]/zgoubi.fai '
      write(*,*) 'Averaging starts at pass # ',j1Pass
      write(*,*) ' ok ? '
      read(*,*)
      
      open(unit=1,file=fileFai)
      read(1,fmt='(a)',end=97,err=98) txt800
      read(1,fmt='(a)',end=97,err=98) txt800
      read(1,fmt='(a)',end=97,err=98) txt800
      read(1,fmt='(a)',end=97,err=98) txt800

      open(unit=2,file='runningAverage_fromFai.out')
      write(2,fmt='(a)') 'i, ipass,s,t, '
     >  //'  sdav , sd2av,'
     >  //'  syav/cm , sy2av,'
     >  //'  stav/mrd , st2av,'
     >  //'  szav/cm , sz2av,'
     >  //'  spav/mrd , sp2av,'
     >  //'  ssav , ss2av,'
     >  //'  stiav/mu_s , sti2av,'
     >  //'  spxav , spx2av,'
     >  //'  spyav , spy2av,'
     >  //'  spzav , spz2av,'
     >  //'  spmav , spm2av, '
     >  //'  sxmi, sxma, '
     >  //'  symi, syma, '
     >  //'  szmi, szma, '
     >  //'  Ekin, Etot, am, Q, G'

      ini =  j1Pass ! pass to start averaging from 
      ifi = 999999
      
      sd = 0.d0
      sy = 0.d0
      st = 0.d0
      sz = 0.d0
      sp = 0.d0
      ss = 0.d0
      st = 0.d0
      sd2 = 0.d0
      sy2 = 0.d0
      st2 = 0.d0
      sz2 = 0.d0
      sp2 = 0.d0
      ss2 = 0.d0
      st2 = 0.d0
      sspx = 0.d0
      sspy = 0.d0
      sspz = 0.d0
      sspm = 0.d0
      sspx2 = 0.d0
      sspy2 = 0.d0
      sspz2 = 0.d0
      sspm2 = 0.d0
      sxmi = 1d10
      sxma = -1d10
      symi = 1d10
      syma = -1d10
      szmi = 1d10
      szma = -1d10
      i = 0
 1        continue
        read(1,fmt='(a)',end=97,err=98) txt800
C        write(*,*) txt800

        read(txt800(646:655),*) ipass
        if(ipass.lt.ini) goto 1
        if(ipass.gt.ifi) goto 10
        read(txt800(125:300),*) d,y,t,z,p,s,ti
        i = i+1
        sd = sd +d ; sd2 = sd2 + d*d
        sy = sy +y ; sy2 = sy2 + y*y
        st = st +t ; st2 = st2 + t*t
        sz = sz +z ; sz2 = sz2 + z*z
        sp = sp +p ; sp2 = sp2 + p*p
        ss = ss +s ; ss2 = ss2 + s*s
        sti = sti +ti ; sti2 = sti2 + ti*ti
        sdav = sd/dble(i) ; sd2av = sd2/dble(i)
        syav = sy/dble(i) ; sy2av = sy2/dble(i)
        stav = st/dble(i) ; st2av = st2/dble(i)
        szav = sz/dble(i) ; sz2av = sz2/dble(i)
        spav = sp/dble(i) ; sp2av = sp2/dble(i)
        ssav = ss/dble(i) ; ss2av = ss2/dble(i)
        stiav = sti/dble(i) ; sti2av = sti2/dble(i)
c        sd = sd +d ; sd2 = sd2 + d*d
c        sy = sy +y ; sy2 = sy2 + y*y
c        st = st +t ; st2 = st2 + t*t

        read(txt800(364:428),*) spx,spy,spz,spmod        
        sspx = sspx +spx ; sspx2 = sspx2 + spx*spx
        sspy = sspy +spy ; sspy2 = sspy2 + spy*spy
        sspz = sspz +spz ; sspz2 = sspz2 + spz*spz
        sspm = sspm +spmod ; sspm2 = sspm2 + spmod*spmod
        spxav = sspx/dble(i) ; spx2av = sspx2/dble(i)
        spyav = sspy/dble(i) ; spy2av = sspy2/dble(i)
        spzav = sspz/dble(i) ; spz2av = sspz2/dble(i)
        spmav = sspm/dble(i) ; spm2av = sspm2/dble(i)

        if(spx .lt. sxmi) sxmi = spx
        if(spx .gt. sxma) sxma = spx
        if(spy .lt. symi) symi = spy
        if(spy .gt. syma) syma = spy
        if(spz .lt. szmi) szmi = spz
        if(spz .gt. szma) szma = spz

        read(txt800(429:544),*)  EK, ET, IT, IRP, SRT, am,Q,G

        s = (c * t)* smu2s
        write(2,fmt='(2(i6,1x),1p,35(e14.6,1x))') i,ipass,s,t,
     >  sdav , sd2av,
     >  syav , sy2av,
     >  stav , st2av,
     >  szav , sz2av,
     >  spav , sp2av,
     >  ssav , ss2av,
     >  stiav , sti2av,
     >  spxav , spx2av, 
     >  spyav , spy2av, 
     >  spzav , spz2av, 
     >  spmav , spm2av,
     >  sxmi, sxma, 
     >  symi, syma, 
     >  szmi, szma,
     >  EK, ET,am,Q,G
      goto 1

 10      continue
      
      goto 99
      
 97      continue
      write(*,*) ' End upon EOF'
      goto 99
 98      continue
      write(*,*) ' End upon read-error'
      
 99      continue
      stop ' End upon read-error'

      end
      
 
