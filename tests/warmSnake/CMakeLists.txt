# Set up warmSnake test

set(test_name "${test_directory}")

configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/${test_directory}/centeredHelix_FIT_save_nofinal_150226.res"
  "${ZGOUBI_TEST_BIN_DIR}/${test_directory}/centeredHelix_FIT_save_nofinal_150226.res"
  COPYONLY
)
configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/${test_directory}/warmSnake.res.expected"
  "${ZGOUBI_TEST_BIN_DIR}/${test_directory}/warmSnake.res.expected"
@ONLY
)

set(ZGOUBI_TEST_IZ 29)
set(ZGOUBI_TEST_ID  3)
set(ZGOUBI_TEST_MMAP 1)
configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/${test_directory}/PARIZ_cmake_configurable.H"
  "${ZGOUBI_TEST_BIN_DIR}/${test_directory}/PARIZ.H"
@ONLY
)
configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/${test_directory}/pariz_cmake_configurable.nml"
  "${ZGOUBI_TEST_BIN_DIR}/${test_directory}/pariz.nml"
@ONLY
)

#set(exts "")
set(res_error 1.0e-09)
set(out_error 1.0e-09)
add_test(
  NAME ${test_name}-test
  COMMAND "${ZGOUBI_TEST_BIN_DIR}/run-test.sh" "${test_name}" "${res_error}" "${out_error}" "rebuild-with-field-map"
)
set_tests_properties(${test_name}-test PROPERTIES PASS_REGULAR_EXPRESSION "Test passed.")
