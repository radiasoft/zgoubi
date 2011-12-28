        subroutine bbkLm(Ptsl,npt,linearmap)
        implicit double precision (a-h,o-z)
        parameter (mxpt=10000)
        dimension Ptsl(9,mxpt)
        double precision  linearmap(6,6) 
        double precision  tmp1,tmp2

        do i = 1, npt
          tmp1 = Ptsl(1,i)
          tmp2 = Ptsl(2,i)
          Ptsl(1,i) = linearmap(1,1)*tmp1+linearmap(1,2)*tmp2
          Ptsl(2,i) = linearmap(2,1)*tmp1+linearmap(2,2)*tmp2
          tmp1 = Ptsl(3,i)
          tmp2 = Ptsl(4,i)
          Ptsl(3,i) = linearmap(3,3)*tmp1+linearmap(3,4)*tmp2
          Ptsl(4,i) = linearmap(4,3)*tmp1+linearmap(4,4)*tmp2
          tmp1 = Ptsl(5,i)
          tmp2 = Ptsl(6,i)
          Ptsl(5,i) = linearmap(5,5)*tmp1+linearmap(5,6)*tmp2
          Ptsl(6,i) = linearmap(6,5)*tmp1+linearmap(6,6)*tmp2                   
        enddo

        return
        end
