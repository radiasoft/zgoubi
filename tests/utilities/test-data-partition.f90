program main
   !! author: Damian Rouson
   !!
   !! summary: verify that data partitioned across images evenly to 
   !! within a difference of one datum between any two images.
   use data_partition_ixfc, only : data_partition
   use assertions_interface, only : assert

   implicit none
   type(data_partition) partition
   integer, parameter :: num_particles=7
   integer my_particles 

   associate( me=>this_image(), ni=>num_images() )
     call partition%define_partition( cardinality=num_particles)
     my_particles = partition%my_last() - partition%my_first() + 1
     associate( quotient=>num_particles/ni, remainder=>mod(num_particles,ni)  )
       call assert( my_particles == quotient + merge(1, 0, me<=remainder), "block distribution" )
     end associate
     call co_sum( my_particles, result_image=1 )
     if (me==1) then
       call assert( my_particles==num_particles, "all particles distributed" )
       print *,"Test passed."
     end if
   end associate

end program
