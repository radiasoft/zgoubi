# Set up test for this directory

set(test_name "${test_directory}")

configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/${test_directory}/${test_name}.res.expected"
  "${ZGOUBI_TEST_BIN_DIR}/${test_directory}/${test_name}.res.expected"
@ONLY
)

configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/${test_directory}/sol30m_761.map"
  "${ZGOUBI_TEST_BIN_DIR}/${test_directory}/sol30m_761.map"
@ONLY
)

set(res_error "1.0e-06")
set(out_error "1.0e-08")
add_test(
  NAME ${test_name}-test
  COMMAND "${ZGOUBI_TEST_BIN_DIR}/run-test.sh" ${test_name} "${res_error}" "${out_error}"
)
set_tests_properties(
  ${test_name}-test
  PROPERTIES PASS_REGULAR_EXPRESSION "Test passed."
)
