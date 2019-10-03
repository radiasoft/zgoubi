# Set up paramerter_namelist reader test

set(test_name "${test_directory}")

add_executable(${test_name} data-partition.f90)
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
