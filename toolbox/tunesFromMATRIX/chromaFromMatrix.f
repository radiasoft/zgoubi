
      implicit double precision (a-h,o-z)

      character(300) tline
      
      open(unit=1,file = 'tunesFromMatrix.out',action='READ')
      open(unit=2,file = 'chromaFromMatrix.out')

      read(1,fmt='(a)') tline
      write(2,*) trim(tline)//'     14    15 '
      read(1,fmt='(a)') tline
      write(2,*) trim(tline)//'  dqx   dqy '
      read(1,fmt='(a)') tline
      write(2,*) trim(tline)
      
      l = 0
      lr = 0
 1    continue
        read(1,fmt='(a)',end=10,err=10) tline
        l = l+1
        read(tline,*) xco, xpco, Qx, Qy, alfx, betx, alfy, bety, Dx,
     >  tof, path, EK, pp0

        lr = lr +1
        if(lr.gt.1) then
          dqx = (qx-qx1)/(pp0-pp1)
          dqy = (qy-qy1)/(pp0-pp1)
          qx1 = qx ; qy1 = qy ; pp1 = pp0
        else
           qx1 = qx ; qy1 = qy ; pp1 = pp0
           dqx = 0.d0 ; dqy = 0.d0
        endif

        write(2,fmt='(1p,15(e15.7,1x))')
     >  xco, xpco, Qx, Qy, alfx, betx, alfy, bety, Dx,
     >  tof, path, EK, pp0, dqx, dqy

      goto 1

 10   continue

      write(*,*) ' Job don. l = ', l, ',   lr = ',lr
      stop
      end
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
