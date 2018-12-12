add_library( interfaces OBJECT
  object_interface.f90
  assertions_interface.F90
  )
add_library( zgoubi_utilities
  assertions_implementation.F90
  $<TARGET_OBJECTS:interfaces>
  )
set_target_properties(zgoubi_utilities
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod")
set_target_properties(interfaces
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod")

target_compile_options(interfaces
  PUBLIC
  "-fcoarray=single")
target_compile_options(zgoubi_utilities
  PUBLIC
  "-fcoarray=single")

if(NO_ASSERTIONS)
  target_compile_definitions( zgoubi_utilities
    PUBLIC "-DUSE_ASSERTIONS=.false."
  )
else()
  target_compile_definitions( zgoubi_utilities
    PUBLIC "-DUSE_ASSERTIONS=.true."
  )
endif()
