# Set up chromPSR test

set(test_name "${test_directory}")

configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/${test_directory}/pariz.nml"
  "${ZGOUBI_TEST_BIN_DIR}/${test_directory}/pariz.nml"
 COPYONLY
)
configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/${test_directory}/${test_name}.res.expected"
  "${ZGOUBI_TEST_BIN_DIR}/${test_directory}/${test_name}.res.expected"
@ONLY
)
configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/${test_directory}/${test_name}.MATRIX.out.expected"
  "${ZGOUBI_TEST_BIN_DIR}/${test_directory}/${test_name}.MATRIX.out.expected"
@ONLY
)

set(exts "MATRIX.out")
set(res_error "4.0e-04")
set(out_error "2.5e-06")
add_test(
  NAME ${test_name}-test
  COMMAND "${ZGOUBI_TEST_BIN_DIR}/run-test.sh" "-v" "-e ${exts}" ${test_name} "${res_error}" "${out_error}"
)
set_tests_properties(
  ${test_name}-test
  PROPERTIES PASS_REGULAR_EXPRESSION "Test passed."
)
