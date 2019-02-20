module numeric_defs
  implicit none

  !private
  public

  ! for absolutely *all* floating-point declarations,
  ! use
  !   real(dbl) :: foo
  ! and similar
  integer, parameter :: dbl = kind(0.d0)  ! double precision

  ! numeric constants
  real(dbl), parameter :: zero_r = 0.0e0_dbl
  real(dbl), parameter :: half =   0.5e0_dbl
  real(dbl), parameter :: one =    1.0e0_dbl
  real(dbl), parameter :: two =    2.0e0_dbl
  real(dbl), parameter :: three =  3.0e0_dbl
  real(dbl), parameter :: four =   4.0e0_dbl
  real(dbl), parameter :: five =   5.0e0_dbl
  real(dbl), parameter :: six =    6.0e0_dbl
  real(dbl), parameter :: seven =  7.0e0_dbl
  real(dbl), parameter :: eight =  8.0e0_dbl
  real(dbl), parameter :: nine =   9.0e0_dbl

  ! simple allocatable arrays -- used for the derived types following
  type iarray1
    integer, allocatable :: r(:)
  end type iarray1
  type rarray1
    real(dbl), allocatable :: r(:)
  end type rarray1

  ! derived types that facilitate non-rectangular arrays
  !     type(array2) nra
  !     allocate(nra%c(ncols))
  !     do i = 1, ncols
  !       allocate(nra%c(i)%r(nrows(i)))
  !     end do
  ! NB: in nra%c(i)%r(j), j is the fast index
  type iarray2
    type(iarray1), allocatable :: c(:)
  end type iarray2
  type rarray2
    type(rarray1), allocatable :: c(:)
  end type rarray2

end module numeric_defs
