C Strats from the longest tracking (mDamp(mxDmp))
C Only takes particles that made it to the end for the max number of damping times.
      implicit double precision (a-h,o-z)
      character(200) txt200, txtF
      character(4) txt4
      parameter (mxDmp = 20)
      dimension mDamp(mxDmp), mPass(mxDmp)
      logical exs, alive, frstFile

      parameter (nBins = 1024)
      dimension alive(nBins)

      data live / 0 /
      data alive, frstFile / nBins*.true., .true. /

C Working hypotheses ----------------------
      data mDamp / 50, 100, 150, 200, 250, 300, 400, 500,
     >12*0/
      data mPass / 
     >22635 ,45270 ,67905 ,90540 ,120721,135810, 181080, 226351,
     >12*0/

      ag1 = 40.3
      ag2 = 40.6
C----------------------------------------       

      open(unit=2,file='polSpectrumFromFai_average.res')
      write(2,*)'# iFile live icount nPass <SX>  <SY>  <SZ>'
     >//'  <|S|> (zgoubi frame)  [a.gamma_min, a.gamma_max]'

      do i = 8, 1, -1
       write(txt4,fmt='(i0)') mDamp(i)
       txtF = 'polSpectrumFromFai_'//trim(txt4)//'damp.out'
       inquire(exist=exs,file=trim(txtF))

         write(*,*) exs, trim(txtF)
            read(*,*)

       if(exs) then

        call system('ln -sf '//trim(txtF)//
     >  ' polSpectrumFromFai_average.in')
        call system('ls -l polSpectrumFromFai_average.in')
          read(*,*)
        nPass = mPass(i)
        open(unit=1,file='polSpectrumFromFai_average.in')
        read(1,fmt='(a)') txt200
        write(*,*) txt200

        kBin = 0
        ssx = 0.d0
        ssy = 0.d0
        ssz = 0.d0
        smod = 0.d0
        mtmp = 0
        icnt = 0
 1      continue
        read(1,fmt='(a)',err=10,end=10) txt200
        read(txt200(1:23),*) ifldr, agam
        read(txt200(37:200),*) pRef, agref,SX,SY,SZ, Smod, mxPass
        mtmp = mxPass

        kBin = kBin + 1
        if(frstFile) then
          alive(kBin) = (mtmp .eq. mPass(i))
          if(alive(kBin)) live = live + 1
        endif

c        write(*,*) trim(txt200)
c        write(*,*)'  kBin: ',mPass(i),mtmp,mPass(i)-mtmp
c        write(*,*)'  kBin: ',kBin,i,alive(kBin)
c             read(*,*)

        if(.not. alive(kBin)) goto 1
        if(agref .lt. ag1) goto 1
        if(agref .gt. ag2) goto 1

        icnt =  icnt +1
           write(*,*) kbin, icnt
        ssx = ssx + SX
        ssy = ssy + SY
        ssz = ssz + SZ
        ssmod = ssmod + Smod
        goto 1

 10     continue
        frstFile = .false.
        close(1)

        write(*,fmt='(a,2(f10.3,1x),a,i0)') 
     >  '# Averages are over '
     >  //'interval [a.gamma_min, a.gamma_max]='
     >  ,ag1,ag2,
     >  ',  at turn number nPass=',nPass
        write(*,*)'# Number of particles survived in first file : ',live
     >  ,' (others are not accounted for in the averaging. '
        write(*,*)'# live nPass <SX>  <SY>  <SZ>  <|S|> (zgoubi frame)'
     >  //'  [a.gamma_min, a.gamma_max]'
        write(*,*) 
     >  live, nPass, ssx/dble(icnt), ssy/dble(icnt), ssz/dble(icnt),
     >  ag1, ag2

        write(2,fmt='(4(i6,1x),1p,3(e14.6,1x),0p,2(f7.2,1x),a)')  
     >  i,live,icnt, nPass, 
     >  ssx/dble(icnt), ssy/dble(icnt), ssz/dble(icnt),
     >  ag1, ag2, trim(txtF)

       endif

      enddo


       close(2)


          stop
          end

