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
      subroutine bbrsp(Ptsl,thispart,a,b,c,sa,sb)
C      subroutine bbrsp(Ptsl,thispart,a,b,c,sa,sb,npt)
      implicit double precision (a-h,o-z)
      parameter (mxpt=10000)
      dimension  Ptsl(9,mxpt) 
      integer  thispart
      double precision  a,b,c,sa,sb
      dimension  rotmat(3,3) 
      double precision  tmp1,tmp2,tmp3

      rotmat(1,1) = 1-sa*(b**2+c**2)
      rotmat(1,2) = a*b*sa+c*sb
      rotmat(1,3) = a*c*sa-b*sb
      rotmat(2,1) = a*b*sa-c*sb
      rotmat(2,2) = 1-sa*(a**2+c**2)
      rotmat(2,3) = b*c*sa+a*sb
      rotmat(3,1) = a*c*sa+b*sb
      rotmat(3,2) = b*c*sa-a*sb
      rotmat(3,3) = 1-sa*(a**2+b**2)

      tmp1 = Ptsl(7,thispart)
      tmp2 = Ptsl(8,thispart)
      tmp3 = Ptsl(9,thispart)

      Ptsl(7,thispart) = 
     > rotmat(1,1)*tmp1+rotmat(1,2)*tmp2+rotmat(1,3)*tmp3
      Ptsl(8,thispart) = 
     > rotmat(2,1)*tmp1+rotmat(2,2)*tmp2+rotmat(2,3)*tmp3
      Ptsl(9,thispart) = 
     > rotmat(3,1)*tmp1+rotmat(3,2)*tmp2+rotmat(3,3)*tmp3
      
      return
      end 
