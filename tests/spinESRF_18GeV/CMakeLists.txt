# Set up spinESRF_18GeV test

configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/spinESRF_18GeV/spinESRF_18GeV.res.expected"
  "${ZGOUBI_TEST_BIN_DIR}/spinESRF_18GeV/spinESRF_18GeV.res.expected"
@ONLY
)
configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/spinESRF_18GeV/spinESRF_18GeV.SPNPRT.Out.expected"
  "${ZGOUBI_TEST_BIN_DIR}/spinESRF_18GeV/spinESRF_18GeV.SPNPRT.Out.expected"
@ONLY
)
set(res_error "1.0e-09") 
set(out_error "1.0e-08") 
add_test(
  NAME spinESRF_18GeV-test 
  COMMAND "${ZGOUBI_TEST_BIN_DIR}/run-test.sh" spinSaturne "${res_error}" "${out_error}"
)
set_tests_properties(
  spinESRF_18GeV-test 
  PROPERTIES PASS_REGULAR_EXPRESSION "Test passed."
)
