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

   associate( me=>this_image(), ni=>num_images() )

     call partition%define_partition( cardinality=num_particles)

     associate( my_first=>partition%my_first(), my_last=>partition%my_last() )

       verify_partitioning: &
       block
         integer my_particles

         my_particles = my_last - my_first + 1

         associate( quotient=>num_particles/ni, remainder=>mod(num_particles,ni)  )
           call assert( my_particles == quotient + merge(1, 0, me<=remainder), "block distribution" )
         end associate

         call co_sum( my_particles, result_image=gatherer )
         if (me==gatherer) call assert( my_particles==num_particles, "all particles distributed" )
       end block verify_partitioning

       verify_gather: &
       block
         real(real64), parameter :: junk=-123456._real64, expected=1._real64
         integer, parameter :: vec_space_dim=3
         real(real64) :: particle_scalar(num_particles), particle_vector(vec_space_dim, num_particles)
         real(real64) :: vector_transpose(num_particles, vec_space_dim)

         particle_scalar(my_first:my_last) = expected !! values to be gathered
         particle_scalar(1:my_first-1)  = junk !! values to be overwritten by the gather
         particle_scalar(my_last+1:)  = junk !! values to be overwritten by the gather

         call partition%gather( particle_scalar, result_image=gatherer )
         select case(me) 
           case(gatherer)
             call assert( all( particle_scalar==expected), "correct values gathered" ) 
               !! expect exact replicas because the only floating point arithmetic involves adding 0.
           case default
             call assert( all( particle_scalar(my_first:my_last)==expected), "source data unchanged" )
         end select

         particle_vector(:, my_first:my_last) = expected !! values to be gathered
         particle_vector(:, 1:my_first-1)  = junk !! values to be overwritten by the gather
         particle_vector(:, my_last+1:)  = junk !! values to be overwritten by the gather

         call partition%gather( particle_vector, result_image=gatherer, dim=2 )

         select case(me) 
           case(gatherer)
             call assert( all( particle_vector==expected), "correct values gathered" ) 
               !! expect exact replicas because the only floating point arithmetic involves adding 0.
           case default
             call assert( all( particle_vector(:, my_first:my_last)==expected), "source data unchanged" )
         end select

         vector_transpose(my_first:my_last, :) = expected !! values to be gathered
         vector_transpose(1:my_first-1, :)  = junk !! values to be overwritten by the gather
         vector_transpose(my_last+1:, :)  = junk !! values to be overwritten by the gather

         call partition%gather( vector_transpose, result_image=gatherer, dim=1 )

         select case(me) 
           case(gatherer)
             call assert( all( vector_transpose==expected), "correct values gathered" ) 
               !! expect exact replicas because the only floating point arithmetic involves adding 0.
           case default
             call assert( all( vector_transpose(my_first:my_last, :)==expected), "source data unchanged" )
         end select

       end block verify_gather

       sync all !! wait for all images' assertions to complete
       if (me==1) print *,"Test passed"
     end associate
   end associate

end program
