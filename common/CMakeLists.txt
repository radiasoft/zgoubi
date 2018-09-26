#-------------------------------------------------------------------------------
#
# Define source files for libzg
# =============================
#
# De-clutter top-level CMakeLists.txt and use modern CMake best practices:
# Treat targets as objects and transitively propagate properties etc.
#
# Provides: ZG_SRCS
#
#-------------------------------------------------------------------------------

# Sources grouped by directory and alphabetized

set(_libzg_srcs
  date2.f
  empty.f
  endjob.f
  gttext.f
  iraz.f
  raz.f
  rffag.f
  rhist.f
  rndm.f
  strcon.f
  strget.f
  time2.f
  vscal.f
  vvect.f
  )
foreach(src IN LISTS _libzg_srcs)
  list(APPEND ZG_SRCS
    "${CMAKE_CURRENT_LIST_DIR}/${src}")
endforeach()
