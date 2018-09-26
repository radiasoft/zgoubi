# Set up spinESRF_18GeV test

configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/spinESRF_18GeV/spinESRF_18GeV.sh.template"
  "${ZGOUBI_TEST_BIN_DIR}/spinESRF_18GeV/spinESRF_18GeV-test.sh"
@ONLY
)
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
add_test(NAME spinESRF_18GeV-test COMMAND "${ZGOUBI_TEST_BIN_DIR}/spinESRF_18GeV/spinESRF_18GeV-test.sh" spinESRF_18GeV)
set_tests_properties(spinESRF_18GeV-test PROPERTIES PASS_REGULAR_EXPRESSION "Test passed.")
