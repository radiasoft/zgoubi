add_library( zgoubi_utilities
  object_interface.f90
  assertions_interface.F90
  assertions_implementation.F90
  data-partition-ixfc.f90
  data-partition-impl.f90
  )
set_target_properties(zgoubi_utilities
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod")
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
