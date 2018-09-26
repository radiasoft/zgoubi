# Set up spinSaturne test

configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/spinSaturne/spinSaturne.sh.template"
  "${ZGOUBI_TEST_BIN_DIR}/spinSaturne/spinSaturne-test.sh"
@ONLY
)
configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/spinSaturne/spinSaturne.res.expected"
  "${ZGOUBI_TEST_BIN_DIR}/spinSaturne/spinSaturne.res.expected"
@ONLY
)
configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/spinSaturne/spinSaturne.SPNPRT.Out.expected"
  "${ZGOUBI_TEST_BIN_DIR}/spinSaturne/spinSaturne.SPNPRT.Out.expected"
@ONLY
)
add_test(NAME spinSaturne-test COMMAND "${ZGOUBI_TEST_BIN_DIR}/spinSaturne/spinSaturne-test.sh" spinSaturne)
set_tests_properties(spinSaturne-test PROPERTIES PASS_REGULAR_EXPRESSION "Test passed.")
