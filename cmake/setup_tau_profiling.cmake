# Find TAU Commander if it's installed, see if the project is setup and set it up if TAU is present but the project is not installed.

# Look for TAU Commander executable on $PATH and in suspected locations
find_program(TAU_COMMANDER
  tau
  HINTS "$ENV{HOME}/taucmdr/bin" "$ENV{HOME}/.local/bin"
  DOC "Path to TAU Commander binary executable")

# Look for .tau project dir, works like .git dirs, usually this should be in ${CMAKE_SOURCE_DIR}
find_path(TAU_PROJECT_DIR
  .tau
  HINTS "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}" "${CMAKE_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}"
  DOC "Location for the TAU project directory")

# Warn user and attempt to install TAU Commander if profiling requested but `tau` not found
if(CMAKE_BUILD_TYPE MATCHES "[Bb]aseline|[Ss]ample|[Cc]ompiler-inst")
  if(NOT TAU_COMMANDER)
    message(WARNING "A profiling build configuration was requested, but TAU Commander was not found. Attempting to install TAU Commander now.")
    execute_process(COMMAND ./install-update-tau.sh
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
    unset(TAU_COMMANDER CACHE)
    find_program(TAU_COMMANDER
      tau
      PATHS "$ENV{HOME}/taucmdr/bin"
      DOC "Path to TAU Commander binary executable")
  endif()
endif()

# If TAU Commander found and profiling requested, but no project directory is present, initialize one
if(TAU_COMMANDER)
  if(NOT TAU_PROJECT_DIR)
    message(STATUS "TAU Commander was found, but no TAU Commander project appears to be initialized")
    if(CMAKE_BUILD_TYPE MATCHES "[Bb]aseline|[Ss]ample|[Cc]ompiler-inst")
      message(WARNING "A profiling build configuration was requested, but no project was initialized. A TAU project will be initialized now.")
      execute_process(COMMAND ./setup_tau_project.sh
	WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
      unset(TAU_PROJECT_DIR CACHE)
      find_path(TAU_PROJECT_DIR
	.tau
	PATHS "${CMAKE_SOURCE_DIR}"
	DOC "Location for the TAU project directory")
      if(NOT TAU_PROJECT_DIR)
	message(FATAL_ERROR "Failed to initialize a TAU Commander project, despite profiling build requested, aborting build.")
      endif()
    endif()
  endif()
endif()

set(TAU_EXE "")
# Make sure that the correct tau measurement is selected
# This may be problematic if multiple targets or applications are present
if(CMAKE_BUILD_TYPE MATCHES "[Bb]aseline|[Ss]ample|[Cc]ompiler-inst")
  string(TOLOWER "${CMAKE_BUILD_TYPE}" TAU_CONFIG)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_CONFIG)
  get_filename_component(TAU_ROOT_DIR "${TAU_PROJECT_DIR}" NAME_WE)
  execute_process(COMMAND "${TAU_COMMANDER}" select ${TAU_CONFIG}
    WORKING_DIRECTORY "${TAU_ROOT_DIR}")

  # Use defaults from the RelWithDebInfo config to get debug symbols and optimization
  set(CMAKE_Fortran_FLAGS_${BUILD_CONFIG} ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}
    CACHE STRING "Fortran compiler flags for the ${CMAKE_BUILD_TYPE} configuration")
  set(CMAKE_EXE_LINKER_FLAGS_${BUILD_CONFIG} ${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO}
    CACHE STRING "Executable linker flags for the ${CMAKE_BUILD_TYPE} configuration")
  set(TAU_EXE "${TAU_COMMANDER} ") # Trailing space is critical! This variable is in a macro expansion
endif()
