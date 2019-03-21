cmake_minimum_required(VERSION 3.10.0)

# Ensure the default configs are defined and select one if none are
# chosen by the user
set ( CMAKE_CONFIGURATION_TYPES "Debug" "RuntimeDebug" "Release" "MinSizeRel" "RelWithDebInfo" )
set ( CMAKE_BUILD_TYPE "Release"
  CACHE STRING "Select which configuration to build." )

# Name project and specify source languages
project(zgoubi
  LANGUAGES Fortran C)

# Print an error message on an attempt to build inside
# the top level of the source tree:
if ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
 message(FATAL_ERROR "ERROR! "
   "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
   " == CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}"
   "\nThis archive does not support in-source builds:\n"
   "To facilitate a correct configuration, please delete the "
   "CMakeCache.txt file and the CMakeFiles/ directory under the "
   "source directory and run something like \n"
   "  $ rm -r CMakeCache.txt CMakeFiles/"
   "\n"
   "Please create a new directory and build in a manner such as"
   "  $ mkdir build\n"
   "  $ cd build\n"
   "  $ FC=gfortran CC=gcc cmake .. -DCMAKE_INSTALL_PREFIX=<install-path>\n"
 )
endif()

set(issues_URL "https://github.com/radiasoft/zgoubi/issues")
if (NOT "${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  message(WARNING "Untested compiler detected: ${CMAKE_Fortran_COMPILER_ID}. Report any issues at ${issues_URL}")
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/maybe_fetch_file.cmake")
find_package(UnixCommands)

# Get map file needed for test
macro(update_verify_warmSnake)
  maybe_fetch_file( tests/warmSnake/warmSnake.map.gz
    "https://sourceforge.net/p/zgoubi/code/HEAD/tree/trunk/exemples/usersGuide/FIT-and-REBELOTE/warmSnake.map.gz?format=raw"
    06fac3c92404100ffcbbaf281136c2e3d565c6091aa24d3a216aae816b012878 # 8 # Real last digit is 8 not 7 use 7 for testing
    have_warm_snake
    )
endmacro()

update_verify_warmSnake()
message( STATUS "have_warm_snake? ${have_warm_snake}")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(common_fflags -g -static-libgfortran -Wall -fbacktrace -pedantic -cpp)
  set(SAVE_EVERYWHERE -fno-automatic)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(common_fflags -g -standard-semantics -fp-model precise -xHOST -traceback -staticlib -cpp)
  set(SAVE_EVERYWHERE -save)
endif()


#======================================================
# Define targets to be built: libraries and executables
#======================================================

find_package(BLAS REQUIRED)
add_library(BLAS UNKNOWN IMPORTED GLOBAL)
list(GET BLAS_LIBRARIES 0 MAIN_BLAS_LIB)
set_target_properties(BLAS
  PROPERTIES IMPORTED_LOCATION "${MAIN_BLAS_LIB}")

message(STATUS "BLAS LIBRARIES: ${BLAS_LIBRARIES}")
message(STATUS "BLAS LINKER_FLAGS: ${BLAS_LINKER_FLAGS}")

list(REMOVE_AT BLAS_LIBRARIES 0)
if(${BLAS_LIBRARIES})
  target_link_libraries(BLAS
    INTERFACE ${BLAS_LIBRARIES})
endif()
if(${BLAS_LINKER_FLAGS})
  target_link_libraries(BLAS
    INTERFACE "${BLAS_LINKER_FLAGS}")
endif()


find_package(LAPACK REQUIRED)
add_library(LAPACK UNKNOWN IMPORTED GLOBAL)
list(GET LAPACK_LIBRARIES 0 MAIN_LAPACK_LIB)
set_target_properties(LAPACK
  PROPERTIES IMPORTED_LOCATION "${MAIN_LAPACK_LIB}")

message(STATUS "LAPACK LIBRARIES: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK LINKER_FLAGS: ${BLAS_LINKER_FLAGS}")

list(REMOVE_AT LAPACK_LIBRARIES 0)
if(${LAPACK_LIBRARIES})
  target_link_libraries(LAPACK
    INTERFACE ${LAPACK_LIBRARIES})
endif()
if(${LAPACK_LINKER_FLAGS})
  target_link_libraries(LAPACK
    INTERFACE "${LAPACK_LINKER_FLAGS}")
endif()

#---------
# libcouplb
#---------
add_library(coupl "")
# Collect & define COUPL_SRCS
include("${CMAKE_CURRENT_SOURCE_DIR}/zgoubi/coupling/CMakeLists.txt")
target_sources(coupl
  PRIVATE ${COUPL_SRCS})
set_source_files_properties(${COUPL_SRCS}
  PROPERTIES COMPILE_OPTIONS ${SAVE_EVERYWHERE})
target_compile_options(coupl PUBLIC ${common_fflags})
target_link_libraries(coupl
  PUBLIC BLAS
  PUBLIC LAPACK)

#------
# libzg
#------
add_library(zg "")
# Collect and define ZG_SRCS, ZG_INCLUDES
include("${CMAKE_CURRENT_SOURCE_DIR}/common/CMakeLists.txt")
target_sources(zg
  PRIVATE ${ZG_SRCS})
set_source_files_properties(${ZG_SRCS}
  PROPERTIES COMPILE_OPTIONS ${SAVE_EVERYWHERE})
target_sources(zg
  PUBLIC ${ZG_INCLUDES})
target_compile_options(zg PUBLIC ${common_fflags})
target_include_directories(zg
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/zgoubi>  # <prefix>/include/zgoubi
)

#------------------
# zgoubi executable
#------------------
add_executable(zgoubi "")
# Collect and define ZGOUBI_SRCS
include("${CMAKE_CURRENT_SOURCE_DIR}/zgoubi/CMakeLists.txt")
target_sources(zgoubi
  PRIVATE ${ZGOUBI_SRCS})
set_source_files_properties(${ZGOUBI_SRCS}
  PROPERTIES COMPILE_OPTIONS ${SAVE_EVERYWHERE})

set( zgoubi_modules
  dynhc.f
  numeric_defs.f90
  taylor.f90
  mathphys_consts.f90
  particle_ixfc.f90
  particle_impl.f90
  dBarrays.f90
)
foreach(module_file ${zgoubi_modules})
  target_sources(zgoubi
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/modules/${module_file}")
endforeach()

target_sources(zgoubi
  PUBLIC ${ZGOUBI_INCLUDES})
target_compile_options(zgoubi PRIVATE ${common_fflags})
target_link_libraries(zgoubi
  PUBLIC zg
  PUBLIC coupl
  PUBLIC zgoubi_utilities)
target_include_directories(zgoubi
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/mod>
  $<INSTALL_INTERFACE:include/zgoubi>  # <prefix>/include/zgoubi
)
set_target_properties( zgoubi
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod"
)

#-------------
# libfminigraf
#-------------
if ( "${CMAKE_SYSTEM_NAME}" MATCHES "Darwin" )
  set(OSNAME DARWIN)
else()
  set(OSNAME LINUX)
  if ( NOT "${CMAKE_SYSTEM_NAME}" MATCHES "Linux" )
    message(WARNING "Untested operating system: ${CMAKE_SYSTEM_NAME}. Report any issues at ${issues_URL}")
  endif()
endif()
add_library(fminigraf "")
target_sources(fminigraf
  PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/zpop/liblns/fminigraf.c")
target_compile_definitions(fminigraf
  PRIVATE -D${OSNAME})
target_compile_options(fminigraf
  PRIVATE -g)

#----------------
# zpop executable
#----------------
add_executable(zpop "")
# Collect and define ZPOP_SRCS
include("${CMAKE_CURRENT_SOURCE_DIR}/zpop/CMakeLists.txt")
target_sources(zpop
  PRIVATE ${ZPOP_SRCS})
set_source_files_properties(${ZPOP_SRCS}
  PROPERTIES COMPILE_OPTIONS ${SAVE_EVERYWHERE})
target_sources(zpop
  PUBLIC ${ZPOP_INCLUDES})
target_compile_options(zpop PRIVATE ${common_fflags})
target_link_libraries(zpop zg fminigraf)
target_include_directories(zpop
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/zgoubi>  # <prefix>/include/zgoubi
)


#==============================================
# Handle installation of binaries and libraries
#==============================================

# Set installation subdirectory for executable files
install(
  TARGETS zgoubi zpop
  DESTINATION bin
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

# Set installation subdirectory for libraries
if(BUILD_SHARED_LIBS)
  install(
    TARGETS zg coupl fminigraf
    DESTINATION lib
    PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
    )
endif()

enable_testing()

set(ZGOUBI_TEST_SRC_DIR ${CMAKE_SOURCE_DIR}/tests )
set(ZGOUBI_TEST_BIN_DIR ${CMAKE_BINARY_DIR}/tests )

include_directories(BEFORE ${ZGOUBI_TEST_BIN_DIR} zgoubi include)

add_subdirectory(utilities)
add_subdirectory(tests)
