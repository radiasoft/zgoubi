submodule(data_partition_ixfc) data_partition_impl
  use assertions_interface, only : assert, assertions
  implicit none

  type(data_partition) singleton
    !! This is the sole instance of this class in the executing image, making 
    !! it the multiverse analogue of the classical Singleton design pattern.

contains

  module procedure define_partitions


    if (allocated(singleton%first_datum)) deallocate(singleton%first_datum)
    if (allocated(singleton%last_datum)) deallocate(singleton%last_datum)

    associate( ni => num_images() )

      call assert( ni<=cardinality, "sufficient data for distribution across images")

      allocate(singleton%first_datum(ni), singleton%last_datum(ni))

      block 
        integer i, image
        do image=1,ni
          associate( remainder => mod(cardinality, ni), quotient => cardinality/ni )
            singleton%first_datum(image) = sum([(quotient+overflow(i, remainder), i=1, image-1)]) + 1
            singleton%last_datum(image) = singleton%first_datum(image) + quotient + overflow(image, remainder) - 1
          end associate
        end do
      end block
    end associate

  contains

    pure function overflow(im, excess) result(extra_datum)
      integer, intent(in) :: im, excess
      integer extra_datum
      extra_datum= merge(1,0,im<=excess)
    end function

  end procedure 

  module procedure first
    if (assertions) call assert( allocated(singleton%first_datum), "allocated(singleton%first_datum)")
    first_datum = singleton%first_datum( image_number )
  end procedure

  module procedure last
    if (assertions) call assert( allocated(singleton%last_datum), "allocated(singleton%last_datum)")
    last_datum = singleton%last_datum( image_number )
  end procedure

  module procedure gather_real_1D_array

    if (present(dim)) call assert (dim==1, "dimensioned partitioned == 1")

    associate( me => this_image() )
      associate( first=>singleton%first(me), last=>singleton%last(me) )
        if (.not. present(result_image)) then
          a(1:first-1)  = 0.
          a(last+1:)  = 0.
          call co_sum(a)
        else
          block
            real(real64), allocatable, dimension(:) :: a_lower, a_upper
            a_lower = a(1:first-1)
            a_upper = a(last+1:)
            a(1:first-1)  = 0.
            a(last+1:)  = 0.
            call co_sum(a, result_image=result_image)
            if (result_image /= me) then
              a(1:first-1) = a_lower
              a(last+1:) = a_upper
            end if
          end block
        end if
      end associate
    end associate
  end procedure

  module procedure gather_real_2D_array

    integer dim_
    if (present(dim)) then
      dim_ = dim
    else
      dim_ = 2
    end if

    associate( me => this_image() )
      associate( first => singleton%first(me), last => singleton%last(me) )
        if (.not. present(result_image)) then
          select case(dim_)
            case(1)
              a(1:first-1, :) = 0.
              a(last+1:, :) = 0.
            case(2)
              a(:, 1:first-1) = 0.
              a(:, last+1:) = 0.
            case default
              error stop "gather_real_2D_array: invalid dim argument"
          end select
          call co_sum(a)
        else
          block
            real(real64), allocatable, dimension(:,:) :: a_lower, a_upper
            select case(dim_)
              case(1)
                a_lower = a(1:first-1, :)
                a_upper = a(last+1:, :)
                a(1:first-1, :) = 0.
                a(last+1:, :) = 0.
              case(2)
                a_lower = a(:, 1:first-1)
                a_upper = a(:, last+1:)
                a(:, 1:first-1) = 0.
                a(:, last+1:) = 0.
              case default
                error stop "gather_real_2D_array: invalid dim argument"
            end select

            call co_sum(a, result_image=result_image)

            if (result_image /= me) then
              select case(dim_)
                case(1)
                  a(1:first-1, :) = a_lower 
                  a(last+1:, :) = a_upper 
                case(2)
                  a(:, 1:first-1) = a_lower 
                  a(:, last+1:) = a_upper 
                case default
                  error stop "gather_real_2D_array: invalid dim argument"
              end select
            end if
          end block
        end if
      end associate
    end associate
  end procedure

end submodule data_partition_impl
