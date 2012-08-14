C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  François Méot
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.
C
C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 51 Franklin Street, Fifth Floor,
C  Boston, MA  02110-1301  USA
C
C  François Méot <fmeot@bnl.gov>
C  Brookhaven National Laboratory     
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
        subroutine normVec(y,num)
        implicit double precision (a-h,o-z)
C S. White & F. Meot, Jan. 2012

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

