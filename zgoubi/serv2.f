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
C  USA
C  -------
      SUBROUTINE SERV2
      use c_ss1_interface, only : JX, JY, NN 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      JY(1)=-2
      JX(1)=-2
      JY(2)=-2
      JX(2)=-1
      JY(3)=-2
      JX(3)=0
      JY(4)=-2
      JX(4)=1
      JY(5)=-2
      JX(5)=2
      JY(6)=-1
      JX(6)=-2
      JY(7)=-1
      JX(7)=-1
      JY(8)=-1
      JX(8)=0
      JY(9)=-1
      JX(9)=1
      JY(10)=-1
      JX(10)=2
      JY(11)=0
      JX(11)=-2
      JY(12)=0
      JX(12)=-1
      JY(13)=0
      JX(13)=0
      JY(14)=0
      JX(14)=1
      JY(15)=0
      JX(15)=2
      JY(16)=1
      JX(16)=-2
      JY(17)=1
      JX(17)=-1
      JY(18)=1
      JX(18)=0
      JY(19)=1
      JX(19)=1
      JY(20)=1
      JX(20)=2
      JY(21)=2
      JX(21)=-2
      JY(22)=2
      JX(22)=-1
      JY(23)=2
      JX(23)=0
      JY(24)=2
      JX(24)=1
      JY(25)=2
      JX(25)=2
      NN=25
      RETURN
      END
