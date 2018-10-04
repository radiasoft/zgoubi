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
set(res_error "1.0e-06")
set(out_error "1.0e-08")
add_test( 
  NAME spinSaturne-test 
  COMMAND "${ZGOUBI_TEST_BIN_DIR}/spinSaturne/spinSaturne-test.sh" spinSaturne "${res_error}" "${out_error}"
)
set_tests_properties(
  spinSaturne-test 
  PROPERTIES PASS_REGULAR_EXPRESSION "Test passed."
)
