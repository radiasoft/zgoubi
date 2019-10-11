submodule(data_partition_ixfc) data_partition_impl
  use assertions_interface, only : assert, assertions
  implicit none

contains

  module procedure define_partition

    integer image

    associate( me => this_image(), ni => num_images() )

      call assert( ni<=cardinality, "sufficient data for distribution across images")

      associate( remainder => mod(cardinality,ni), quotient => cardinality/ni )
        singleton%my_first_datum = sum([(quotient+overflow(image,remainder), image=1, me-1)]) + 1
        singleton%my_last_datum = singleton%my_first_datum + quotient + overflow(me,remainder) - 1
      end associate
    end associate

  contains

    pure function overflow(image, excess) result(extra_datum)
      integer, intent(in) :: image, excess
      integer extra_datum
      extra_datum= merge(1,0,image<=excess)
    end function

  end procedure 

  module procedure my_first
   if (assertions) call assert( allocated(singleton%my_first_datum), "partition lower bound defined")
   first_datum = singleton%my_first_datum
  end procedure

  module procedure my_last
   if (assertions) call assert( allocated(singleton%my_last_datum), "partition upper bound defined")
   last_datum = singleton%my_last_datum
  end procedure

  module procedure gather_real_1D_array
    associate( my_first=>singleton%my_first(), my_last=>singleton%my_last() )
      a(1:my_first-1)  = 0.
      a(my_last+1:)  = 0.
      call co_sum(a,result_image)
    end associate
  end procedure

  module procedure gather_real_2D_array
  end procedure

end submodule data_partition_impl
