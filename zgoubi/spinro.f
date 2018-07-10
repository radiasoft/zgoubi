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
C  Francois Meot <fmeot@bnl.gov>
C  Brookhaven National Laboratory  
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE SPINRO(ANGLE,V,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension x(4), v(3)
      parameter(pi=4.d0*atan(1.d0))
      parameter(dtr=pi/180.d0)
      V0=dsqrt(V(1)**2 + V(2)**2 + V(3)**2)
      angle_tmp=V0*ANGLE/2.d0
      if(angle_tmp .ge. 360.d0) then
         angle_tmp=angle-INT(angle_tmp/360.d0)*360d0
      endif
      if(V0 .gt. 0.d0) then
         X(1)=cos(dtr*angle_tmp)
         SA=sin(dtr*angle_tmp)/V0
         X(2)=V(1)*SA
         X(3)=V(2)*SA
         X(4)=V(3)*SA
      else
         X(1)=1.d0
         X(2)=0.d0
         X(3)=0.d0
         X(4)=0.d0
      endif
c      write(101,*) angle, V0, angle_tmp,X(1),X(2),X(3),X(4)
      return
      end
