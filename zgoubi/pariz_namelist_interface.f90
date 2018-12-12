module pariz_namelist_interface
   !! author: Damian Rouson
   !!
   !! Export global variables and a namelist file reader
   implicit none

   private
   public :: IZ, ID, MMAP, MXX, MXY
   public :: initialize_input_parameters
   public :: all_namelist_values_read

   integer, parameter :: unread=-1
     !! invalid value for distinguishing default-initialized namelist variables
     !! from namelist variables that have been read from the corresponding namelist

   integer, protected :: IZ=unread, ID=unread, MMAP=unread, MXX=unread, MXY=unread

   interface

     module subroutine initialize_input_parameters( pariz_namelist_file )
       !! Initialize the global variables with values read from the pariz namelist
       !! in a 'pariz.nml' file
       implicit none
       character(len=*), intent(in) :: pariz_namelist_file
     end subroutine

     pure module subroutine all_namelist_values_read
     end subroutine

   end interface

end module
