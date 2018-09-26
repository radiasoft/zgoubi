# Set up warmsnake test

configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/warmSnake/centeredHelix_FIT_save_nofinal_150226.res"
  "${ZGOUBI_TEST_BIN_DIR}/warmSnake/centeredHelix_FIT_save_nofinal_150226.res"
  COPYONLY
)
configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/warmSnake/zgoubi-integration-test.sh.template"
  "${ZGOUBI_TEST_BIN_DIR}/warmSnake/zgoubi-integration-test.sh"
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

add_test(NAME zgoubi-integration COMMAND "${ZGOUBI_TEST_BIN_DIR}/warmSnake/zgoubi-integration-test.sh")
set_tests_properties(zgoubi-integration PROPERTIES PASS_REGULAR_EXPRESSION "Test passed.")
