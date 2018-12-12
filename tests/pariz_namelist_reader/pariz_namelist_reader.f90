program main
   !! author: Damian Rouson
   !!
   !! Verify that the pariz_namelist reader completes
   use pariz_namelist_interface, only : initialize_input_parameters
   implicit none

   call initialize_input_parameters('pariz.nml') !! open and read a pariz_namelist

   print *,"Test passed."

end program
