submodule(c_ss1_interface) c_ss1_implementation
  implicit none

contains

  module procedure ensure_xyz_allocation
     call allocate_if_not_allocated(X, MXX)
     call allocate_if_not_allocated(Y, MXY)
     call allocate_if_not_allocated(Z, IZ)
  contains
    subroutine allocate_if_not_allocated(array,array_size)
      use assertions_interface, only : assert, assertions
      real(real64), intent(inout), allocatable :: array(:)
      integer, intent(in) :: array_size
      integer, parameter :: success=0
      integer istat
      if (.not. allocated(X)) then
        allocate(array(array_size), stat=istat)
        if (assertions) then
          call assert(istat==success, "ensure_xyz_allocation: array allocation failed")
        end if
      end if 
    end subroutine
  end procedure

end submodule
