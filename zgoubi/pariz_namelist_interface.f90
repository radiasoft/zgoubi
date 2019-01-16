module pariz_namelist_interface
   !! author: Damian Rouson
   !!
   !! Export global variables and a namelist file reader
   implicit none

   private
   public :: IZ, ID, MMAP, MXX, MXY
   public :: initialize_input_parameters

   integer, protected :: IZ=29, ID=3, MMAP=8, MXX=801, MXY=29
     !! set default values for the input parameters

   interface

     module subroutine initialize_input_parameters( pariz_namelist_file )
       !! Initialize the global variables with values read from the pariz namelist
       !! in a 'pariz.nml' file
       implicit none
       character(len=*), intent(in) :: pariz_namelist_file
     end subroutine

   end interface

end module
