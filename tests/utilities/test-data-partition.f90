program main
   !! author: Damian Rouson
   !!
   !! summary: verify data partitioning across images and data gathering
   use data_partition_ixfc, only : data_partition
   use assertions_interface, only : assert
   use iso_fortran_env, only : real64

   implicit none
   type(data_partition) partition
   integer, parameter :: num_particles=31, gatherer=1, num_steps=9
   integer iteration

   call partition%define_partitions( cardinality=num_particles)

   associate( me=>this_image() )

     associate( my_first=>partition%first(me), my_last=>partition%last(me) )
       do iteration=1,num_steps
         call verify_partitioning( me, my_first, my_last )
         call verify_all_gather( my_first, my_last )
         call verify_gather( me, my_first, my_last )
       end do
     end associate

     sync all !! wait for all images' assertions to complete
     if (me==1) print *,"Test passed"

   end associate

contains 

  subroutine output_2D_array(prefix, array) ! output for debugging
    character(len=*), intent(in) :: prefix
    real(real64), intent(in) :: array(:,:)
    integer m
    do m=1,num_particles
      print *, prefix // " gather, image ",this_image(),", particle_vector(:, ",m,"): ", array(:,m)
    end do 
  end subroutine

  subroutine verify_partitioning( image, first, last )
    !! Verify that the data is partitioned across images evenly to 
    !! within a difference of one datum between any two images.
    integer, intent(in) :: first, last, image
    integer particles
    particles = last - first + 1

    associate( ni=>num_images() )
      associate( quotient=>num_particles/ni, remainder=>mod(num_particles,ni)  )
        call assert( particles == quotient + merge(1, 0, image<=remainder), "block distribution" )
      end associate
    end associate

    call co_sum(particles)
    call assert( particles==num_particles, "all particles distributed" )
  end subroutine

  subroutine verify_all_gather(first, last) 
    integer, intent(in) :: first, last
    real(real64), parameter :: junk=-12345._real64, expected=1._real64
    integer, parameter :: vec_space_dim=3
    real(real64) :: particle_scalar(num_particles), particle_vector(vec_space_dim, num_particles)
    real(real64) :: vector_transpose(num_particles, vec_space_dim)

    particle_scalar(first:last) = expected !! values to be gathered
    particle_scalar(1:first-1)  = junk !! values to be overwritten by the gather
    particle_scalar(last+1:)  = junk !! values to be overwritten by the gather

    call partition%gather(particle_scalar)
    
    call assert( all( particle_scalar==expected), "all( particle_scalar==expected)" ) 

    particle_vector(:, first:last) = expected !! values to be gathered
    particle_vector(:, 1:first-1)  = junk !! values to be overwritten by the gather
    particle_vector(:, last+1:)  = junk !! values to be overwritten by the gather

    call partition%gather(particle_vector)

    call assert( all( particle_vector==expected), "all( particle_vector==expected)" ) 

    vector_transpose(first:last, :) = expected !! values to be gathered
    vector_transpose(1:first-1, :)  = junk !! values to be overwritten by the gather
    vector_transpose(last+1:, :)  = junk !! values to be overwritten by the gather

    call partition%gather( vector_transpose, dim=1)

    call assert( all( vector_transpose==expected), "all( vector_transpose==expected)" ) 
  end subroutine verify_all_gather

  subroutine verify_gather( image, first, last )
    integer, intent(in) :: image, first, last
    real(real64), parameter :: junk=-12345._real64, expected=1._real64
    integer, parameter :: vec_space_dim=3, gatherer=1
    real(real64) :: particle_scalar(num_particles), particle_vector(vec_space_dim, num_particles)
    real(real64) :: vector_transpose(num_particles, vec_space_dim)

    particle_scalar(first:last) = expected !! values to be gathered
    particle_scalar(1:first-1)  = junk !! values to be preseved except on the gather's result_image
    particle_scalar(last+1:)  = junk !! values to be preserved except on the gather's result_image

    call partition%gather(particle_scalar, result_image=gatherer)
    
    if (image==gatherer) then
      call assert( all( particle_scalar==expected), "all( particle_scalar==expected)" )
    else
      call assert( all( particle_scalar(1:first-1)==junk), "verify_gather: lower scalar data unchanged)" )
      call assert( all( particle_scalar(first:last)==expected), "verify_gather: expected scalar data gathered" )
      call assert( all( particle_scalar(last+1:)==junk), "verify_gather: upper scalar data unchanged)" )
    end if

    particle_vector(:, first:last) = expected !! values to be gathered
    particle_vector(:, 1:first-1)  = junk !! values to be overwritten by the gather
    particle_vector(:, last+1:)  = junk !! values to be overwritten by the gather

    call partition%gather(particle_vector, result_image=gatherer)

    if (image==gatherer) then
      call assert( all( particle_vector==expected), "all( particle_vector==expected)" )
    else
      call assert( all( particle_vector(:,1:first-1)==junk), "verify_gather: lower vector data unchanged)")
      call assert( all( particle_vector(:,first:last)==expected), "verify_gather: expected vector data gathered" )
      call assert( all( particle_vector(:,last+1:)==junk), "verify_gather: upper vector data unchanged)" )
    end if

    vector_transpose(first:last, :) = expected !! values to be gathered
    vector_transpose(1:first-1, :)  = junk !! values to be overwritten by the gather
    vector_transpose(last+1:, :)  = junk !! values to be overwritten by the gather

    call partition%gather( vector_transpose, result_image=1, dim=1)

    if (image==gatherer) then
      call assert( all( vector_transpose==expected), "all( particle_vector==expected)" ) 
    else
      call assert( all( vector_transpose(1:first-1,:)==junk), "verify_gather: lower transpose data unchanged)" )
      call assert( all( vector_transpose(first:last,:)==expected), "verify_gather: expected transpose data gathered" )
      call assert( all( vector_transpose(last+1:,:)==junk), "verify_gather: upper transpose data unchanged)" )
    end if

  end subroutine verify_gather

end program
