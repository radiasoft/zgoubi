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

include_guard(GLOBAL)
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

set(_libzg_includes
  MXLD.H
  C.DON.H
  C.CDF.H)
foreach(include IN LISTS _libzg_includes)
  # Note: This should be `include()`-ed from the top level Zgoubi CMakeLists.txt file
  # ${CMAKE_CURRENT_SOURCE_DIR} will evaluate to the top level Zgoubi dir in this case
  list(APPEND ZG_INCLUDES
    "${CMAKE_CURRENT_SOURCE_DIR}/include/${include}")
endforeach()
