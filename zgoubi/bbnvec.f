        subroutine normVec(y,num)
        implicit double precision (a-h,o-z)

        parameter (mxpt=10000)
        dimension y(2,mxpt)
        double precision  twopi,epsilon
        dimension  x1(mxpt),x2(mxpt)

        epsilon = 1.0d-18

        twopi = 4.0d0*dasin(1.0d0)
        call random_number(x2)
        call random_number(x1)

        do i = 1, num
          if(x1(i).eq.0.0d0) x1(i) = epsilon
          y(1,i) = sqrt(-2.0*log(x1(i)))*dcos(twopi*x2(i))
          y(2,i) = sqrt(-2.0*log(x1(i)))*dsin(twopi*x2(i))
c         write(*,*) '---------------------'
c         write(*,*) i, y(1,i),y(2,i)
c         write(*,*) '---------------------'
        enddo
        return
        end

