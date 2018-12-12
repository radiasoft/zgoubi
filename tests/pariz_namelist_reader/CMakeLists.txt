# Set up paramerter_namelist reader test

# Copy input file to build tree
configure_file(
  "${CMAKE_SOURCE_DIR}/zgoubi/pariz.nml"
  "${CMAKE_CURRENT_BINARY_DIR}/pariz.nml"
  COPYONLY
)

set(test_name "${test_directory}")


add_executable( ${test_name}
  pariz_namelist_reader.f90
  ${CMAKE_SOURCE_DIR}/zgoubi/pariz_namelist_interface.f90
  ${CMAKE_SOURCE_DIR}/zgoubi/pariz_namelist_implementation.f90
)
target_link_libraries(${test_name} zgoubi_utilities)

set_target_properties( ${test_name}
   PROPERTIES
   Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod"
)

add_test(
  NAME ${test_name}-test
  COMMAND "${test_name}"
)
set_tests_properties(${test_name}-test PROPERTIES PASS_REGULAR_EXPRESSION "Test passed.")
