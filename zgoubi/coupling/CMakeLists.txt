#-------------------------------------------------------------------------------
#
# Define source files for libcoupl
# ================================
#
# De-clutter top-level CMakeLists.txt and use modern CMake best practices:
# Treat targets as objects and transitively propagate properties etc.
#
# Provides: COUPL_SRCS
#
#-------------------------------------------------------------------------------

# Sources grouped by directory and alphabetized
set( _libcoupl_srcs
  tunesc/extrac.f
  tunesc/hamilt.f
  tunesc/positi.f
  tunesc/propag.f
  tunesc/propag.f

  libraries/libeis/balanc.f
  libraries/libeis/balbak.f
  libraries/libeis/cdiv.f
  libraries/libeis/elmhes.f
  libraries/libeis/eltran.f
  libraries/libeis/hqr.f
  libraries/libeis/hqr2.f
  libraries/libeis/rg.f
  libraries/libeis/rg.f

  libraries/libdet/getdet/getdet.f

  libraries/libinv/determinant-1/lin/cgedi.f
  )

foreach(src IN LISTS _libcoupl_srcs)
  list(APPEND COUPL_SRCS
    "${CMAKE_CURRENT_LIST_DIR}/${src}")
endforeach()
