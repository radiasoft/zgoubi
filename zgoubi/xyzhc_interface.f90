!  ZGOUBI, a program for computing the trajectories of charged particles
!  in electric and magnetic fields
!  Copyright (C) 1988-2007  François Mot
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor,
!  Boston, MA  02110-1301  USA
!
!  François Mot <fmeot@bnl.gov>
!  BNL
!  C-AD, Bldg 911
!  Upton, NY, 11973
!  USA
!  -------
module xyzhc_interface
  !! author: Damian Rouson
  !! 
  !! Export common global array and scalar variables
  use iso_fortran_env, only : real64
  implicit none

  private
  public :: XH, YH, ZH
  public :: IXMA,JYMA,KZMA 
  public :: ensure_xyzhc_allocation

  real(real64), allocatable, dimension(:), save :: XH, YH, ZH
  integer, save :: IXMA, JYMA, KZMA 

  interface

    module subroutine ensure_xyzhc_allocation(MXX,MXY,IZ)
      !! Allocate XH, YH, and ZH if and only if they are currently unallocated
      implicit none
      integer, intent(in) :: MXX, MXY, IZ
    end subroutine
  
  end interface
   
end module
