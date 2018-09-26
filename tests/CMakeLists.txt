set( test_directories
  warmSnake
  spinESRF_18GeV
  spinSaturne
)
foreach(test_directory ${test_directories})
  add_subdirectory(${test_directory})
endforeach()
