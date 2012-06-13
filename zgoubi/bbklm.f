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
