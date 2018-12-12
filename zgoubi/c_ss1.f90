module c_ss1
  use iso_fortran_env, only : real64
  implicit none

  real(real64), allocatable, save :: X(:), Y(:), Z(:)
  integer JY,JX,JZ,NX,NY,NZ,NN

  COMMON/SS1/JY(25),JX(25),JZ(25),NX,NY,NZ,NN

contains

  subroutine allocate_X_Y_if_unallocated(MXX,MXY)
     integer, intent(in) :: MXX,MXY
     if (.not. allocated(X)) allocate(X(MXX))
     if (.not. allocated(Y)) allocate(Y(MXY))
  end subroutine

  subroutine allocate_X_Y_Z_if_unallocated(MXX,MXY,IZ)
     integer, intent(in) :: MXX,MXY,IZ
     if (.not. allocated(X)) allocate(X(MXX))
     if (.not. allocated(Y)) allocate(Y(MXY))
     if (.not. allocated(Z)) allocate(X(IZ))
  end subroutine

end module
