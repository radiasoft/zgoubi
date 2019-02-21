module numeric_defs
  implicit none

  private
  public :: dbl
  public :: iarray1, iarray2, rarray1, rarray2
   !! ragged-edged array types
  public :: zero_r, half, one, two, three, four, five, six, seven, eight, nine
   !! numerical constants

  integer, parameter :: dbl = kind(0.d0)  !! machine-specific double precision
    !! for absolutely *all* floating-point declarations,
    !! use
    !!   real(dbl) :: foo
    !! and similar

  ! Numeric constants
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

  type iarray1
    !! encapsulated integer allocatabl array for creating ragged-edged arrays arrays
    integer, allocatable :: r(:)
  end type iarray1
  type rarray1
    !! encapsulated real allocatabl array for creating ragged-edged arrays arrays
    real(dbl), allocatable :: r(:)
  end type rarray1

  type iarray2
    !! Usage:
    !!     type(array2) nra
    !!     allocate(nra%c(ncols))
    !!     do i = 1, ncols
    !!       allocate(nra%c(i)%r(nrows(i)))
    !!     end do
    !! NB: in nra%c(i)%r(j), j is the fast index
    type(iarray1), allocatable :: c(:)
  end type iarray2
  type rarray2
    type(rarray1), allocatable :: c(:)
  end type rarray2

end module numeric_defs
