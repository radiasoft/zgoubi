
      pi = 4. * atan(1.)

      N = 200
      do i = 1, N+1
        f1 = float(i-1) /float(N) * 2. * pi
        do j = 1, N+1
          f2 = float(j-1) / float(N) * 2. * pi
          if(abs(cosh(f1)*cos(f2)) .le. 1. 
     >    .and.
     >    abs(cos(f1)*cosh(f2)) .le. 1. ) then
            write(88,*) f1/(2.*pi),f2/(2.*pi),' 1' 
          else
            write(88,*) f1/(2.*pi),f2/(2.*pi),' 0'
          endif       
        enddo       
      enddo
              stop
               end
