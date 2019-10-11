program main
   !! author: Damian Rouson
   !!
   !! summary: verify that data partitioned across images evenly to 
   !! within a difference of one datum between any two images.
   use data_partition_ixfc, only : data_partition
   use assertions_interface, only : assert
   use iso_fortran_env, only : real64

   implicit none
   type(data_partition) partition
   integer, parameter :: num_particles=7, gatherer=1
   integer my_particles, i 
   real(real64), parameter :: junk=-123456._real64, expected=1._real64
   real(real64) :: particle_data(num_particles)

   associate( me=>this_image(), ni=>num_images() )
     call partition%define_partition( cardinality=num_particles)

     associate( my_first=>partition%my_first(), my_last=>partition%my_last() )
       my_particles = my_last - my_first + 1

       associate( quotient=>num_particles/ni, remainder=>mod(num_particles,ni)  )
         call assert( my_particles == quotient + merge(1, 0, me<=remainder), "block distribution" )
       end associate

       call co_sum( my_particles, result_image=1 )
       if (me==gatherer) call assert( my_particles==num_particles, "all particles distributed" )

       particle_data(my_first:my_last) = expected !! values to be gathered
       particle_data(1:my_first-1)  = junk !! values to be overwritten by the gather
       particle_data(my_last+1:)  = junk !! values to be overwritten by the gather

       call partition%gather( particle_data, result_image=1 )
       select case(me) 
         case(gatherer)
           call assert( all( particle_data==expected), "correct values gathered" ) 
             !! exact replicas expected because there is no floating point arithmetic involved
         case default
           call assert( all( particle_data(my_first:my_last)==expected), "source data unchanged" )
       end select
       sync all !! wait for all images' assertions to complete
       if (me==1) print *,"Test passed"
     end associate
   end associate
  

end program
