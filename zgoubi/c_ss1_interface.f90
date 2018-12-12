module c_ss1_interface
  use iso_fortran_env, only : real64
  implicit none

  private
  public :: X, Y, Z
  public :: JY, JX, JZ, NX, NY, NZ, NN
  public :: ensure_xyz_allocation

  real(real64), allocatable, save :: X(:), Y(:), Z(:)
  integer JY(25), JX(25), JZ(25), NX, NY, NZ, NN

  interface

    module subroutine ensure_xyz_allocation(MXX, MXY, IZ)
      !! allocate X, Y, and Z if not already allocated
      implicit none
      integer, intent(in) :: MXX, MXY, IZ
    end subroutine

  end interface

end module
