module data_partition_ixfc
  !! distribute data identification numbers across images such that the number of 
  !! items differs by at most 1 between any two images.
  use iso_fortran_env, only : real64
  implicit none

  private
  public :: data_partition

  type data_partition
    !! encapsulate a description of the data subset the executing image owns
    private
    integer, allocatable :: first_datum(:), last_datum(:)
  contains
    procedure, nopass :: define_partitions
    procedure, nopass :: first
    procedure, nopass :: last
    procedure, nopass, private :: gather_real_2D_array, gather_real_1D_array
    generic :: gather => gather_real_2D_array, gather_real_1D_array
  end type

  interface 

    module subroutine define_partitions(cardinality)
      !! define the range of data identification numbers owned by the executing image
      integer, intent(in) :: cardinality
    end subroutine

    pure module function first(image_number) result(first_datum)
      !! the result is the first identification number owned by the executing image
      implicit none
      integer, intent(in) :: image_number
      integer first_datum
    end function

    pure module function last(image_number) result(last_datum)
      !! the result is the last identification number owned by the executing image
      implicit none
      integer, intent(in) :: image_number
      integer last_datum
    end function

    !! Gathers are inherently expensive and are best used either
    !! 1. Near the beginning/end of execution to amortize costs across an entire run or
    !! 2. Temporarily while developing/debugging code.

    module subroutine gather_real_1D_array( a, result_image, dim )
      !! Gather the elements of an 1D array distributed along dimension dim onto result_image
      real(real64), intent(inout) :: a(:)
      integer, intent(in), optional :: result_image
      integer, intent(in), optional :: dim
    end subroutine

    module subroutine gather_real_2D_array( a, result_image, dim )
      !! Gather the elements of an 2D array distributed along dimension dim onto result_image
      real(real64), intent(inout) :: a(:,:)
      integer, intent(in), optional :: result_image
      integer, intent(in), optional :: dim
    end subroutine

  end interface

end module data_partition_ixfc
