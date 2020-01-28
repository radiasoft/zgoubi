configure_file(
  "${ZGOUBI_TEST_SRC_DIR}/run-test.sh.template"
  "${ZGOUBI_TEST_BIN_DIR}/run-test.sh"
@ONLY
)

set( test_directories
  warmSnake
  chromPSR
  spinESRF_18GeV
  spinSaturne
  solenoid
  pariz_namelist_reader
  FITandTWISS
  SVDOC_Saturne
)

foreach(test_directory ${test_directories})
  add_subdirectory(${test_directory})
endforeach()
