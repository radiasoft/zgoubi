module allocations_interface
  implicit none

  private
  public :: allocate_if_unallocated

  interface allocate_if_unallocated
    module procedure allocate_1D_character_array
    module procedure allocate_2D_character_array
    module procedure allocate_1D_real_array
    module procedure allocate_2D_real_array
    module procedure allocate_1D_integer_array
  end interface

  interface

    module subroutine allocate_1D_character_array(array,array_size)
      implicit none
      character(len=*), allocatable, intent(inout):: array(:)
      integer, intent(in) :: array_size
    end subroutine

    module subroutine allocate_2D_character_array(array,rows,columns)
      implicit none
      character(len=*), allocatable, intent(inout):: array(:,:)
      integer, intent(in) :: rows, columns
    end subroutine

    module subroutine allocate_1D_real_array(array,array_size)
      use iso_fortran_env, only : real64
      implicit none
      real(real64), allocatable, intent(inout):: array(:)
      integer, intent(in) :: array_size
    end subroutine

    module subroutine allocate_2D_real_array(array,rows,columns)
      use iso_fortran_env, only : real64
      implicit none
      real(real64), allocatable, intent(inout):: array(:,:)
      integer, intent(in) :: rows, columns
    end subroutine

    module subroutine allocate_1D_integer_array(array,array_size)
      implicit none
      integer, allocatable, intent(inout):: array(:)
      integer, intent(in) :: array_size
    end subroutine

    module subroutine allocate_2D_integer_array(array,rows,columns)
      implicit none
      integer, allocatable, intent(inout):: array(:,:)
      integer, intent(in) :: rows, columns
    end subroutine

  end interface

end module
