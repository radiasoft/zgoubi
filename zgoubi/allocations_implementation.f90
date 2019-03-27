submodule(allocations_interface) allocations_implementation
  use assertions_interface, only : assert, assertions
  implicit none

contains

    module procedure allocate_1D_character_array
      integer istat
      if (.not. allocated(array)) then 
        allocate(array(array_size), stat=istat)
        if (assertions) call assert(istat==0,"array allocation succeeded")
      end if 
    end procedure

    module procedure allocate_2D_character_array
      integer istat
      if (.not. allocated(array)) then 
        allocate(array(rows,columns), stat=istat)
        if (assertions) call assert(istat==0,"array allocation succeeded")
      end if 
    end procedure

    module procedure allocate_1D_real_array
      integer istat
      if (.not. allocated(array)) then 
        allocate(array(array_size), stat=istat)
        if (assertions) call assert(istat==0,"array allocation succeeded")
      end if 
    end procedure

    module procedure allocate_2D_real_array
      integer istat
      if (.not. allocated(array)) then 
        allocate(array(rows,columns), stat=istat)
        if (assertions) call assert(istat==0,"array allocation succeeded")
      end if 
    end procedure

    module procedure allocate_1D_integer_array
      integer istat
      if (.not. allocated(array)) then 
        allocate(array(array_size), stat=istat)
        if (assertions) call assert(istat==0,"array allocation succeeded")
      end if 
    end procedure

    module procedure allocate_2D_integer_array
      integer istat
      if (.not. allocated(array)) then 
        allocate(array(rows,columns), stat=istat)
        if (assertions) call assert(istat==0,"array allocation succeeded")
      end if 
    end procedure

end submodule
