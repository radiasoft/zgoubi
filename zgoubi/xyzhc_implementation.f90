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
submodule(xyzhc_interface) xyzhc_implementation
  implicit none
contains

  module procedure ensure_xyzhc_allocation
    use assertions_interface, only : assert, assertions
    integer allocation_status
    integer, parameter :: success=0

    if (.not. allocated(XH)) allocate(XH(MXX),stat=allocation_status)
    if (assertions) call assert(allocation_status==success,"ensure_xyzhc_allocation: XH allocation succeeded")

    if (.not. allocated(YH)) allocate(YH(MXY),stat=allocation_status)
    if (assertions) call assert(allocation_status==success,"ensure_xyzhc_allocation: YH allocation succeeded")

    if (.not. allocated(ZH)) allocate(ZH(IZ),stat=allocation_status)
    if (assertions) call assert(allocation_status==success,"ensure_xyzhc_allocation: ZH allocation succeeded")
  end procedure
   
end submodule
