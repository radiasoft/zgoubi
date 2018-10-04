# Set up warmsnake test

configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/warmSnake/centeredHelix_FIT_save_nofinal_150226.res"
  "${ZGOUBI_TEST_BIN_DIR}/warmSnake/centeredHelix_FIT_save_nofinal_150226.res"
  COPYONLY
)
configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/warmSnake/warmSnake-test.sh.template"
  "${ZGOUBI_TEST_BIN_DIR}/warmSnake/warmSnake-test.sh"
@ONLY
)
configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/warmSnake/zgoubi.res.expected"
  "${ZGOUBI_TEST_BIN_DIR}/warmSnake/zgoubi.res.expected"
@ONLY
)

set(ZGOUBI_TEST_IZ 29)
set(ZGOUBI_TEST_ID  3)
set(ZGOUBI_TEST_MMAP 1)
configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/warmSnake/PARIZ_cmake_configurable.H"
  "${ZGOUBI_TEST_BIN_DIR}/warmSnake/PARIZ.H"
@ONLY
)

set(res_error 1.0e-09)
set(out_error 1.0e-09)
add_test(
  NAME warmSnake 
  COMMAND "${ZGOUBI_TEST_BIN_DIR}/warmSnake/warmSnake-test.sh" warmSnake "${res_error}" "${out_error}"
)
set_tests_properties(warmSnake PROPERTIES PASS_REGULAR_EXPRESSION "Test passed.")
