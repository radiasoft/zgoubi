set(test "data-partition")
set(executable "test-${test}")
set(source "${executable}.f90" )

add_executable(${executable} ${source})
target_link_libraries(${executable} zgoubi_utilities)
set_target_properties(${executable}
   PROPERTIES
   Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod"
)

set(num_test_images 3)
add_caf_test( ${test} ${num_test_images} ${CMAKE_CURRENT_BINARY_DIR}/${executable})
  
  
  #WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
