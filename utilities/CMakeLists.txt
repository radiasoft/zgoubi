add_library( interfaces OBJECT
  object_interface.f90
  assertions_interface.F90
  data-partition-impl.f90
  data-partition-ixfc.f90
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
target_link_libraries(interfaces
  PUBLIC OpenCoarrays::caf_mpi_static)
target_link_libraries(zgoubi_utilities
  PUBLIC OpenCoarrays::caf_mpi_static)
if(NO_ASSERTIONS)
  target_compile_definitions( zgoubi_utilities
    PUBLIC "-DUSE_ASSERTIONS=.false."
  )
else()
  target_compile_definitions( zgoubi_utilities
    PUBLIC "-DUSE_ASSERTIONS=.true."
  )
endif()
