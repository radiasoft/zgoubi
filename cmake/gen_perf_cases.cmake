# This file is intended to be invoked from CMake in script mode.
# This is used to emulate the envsubst command in gen_template.sh
# found at:
#     https://github.com/radiasoft/rszgoubi/tree/master/examples/SRDampingInESRFRing/coupled_template

if(NOT DEFINED COUPLING_TEMPLATE)
  set( COUPLING_TEMPLATE coupling0.58.template )
endif()
if(NOT DEFINED OUT_FILE_BASENAME)
  set( OUT_FILE_BASENAME coupling )
endif()
if(NOT DEFINED NUM_RUNS)
  set( NUM_RUNS 11)
endif()
if(NOT DEFINED NUM_PARTICLES)
  set( NUM_PARTICLES 1100)
endif()

configure_file(${COUPLING_TEMPLATE} ${OUT_FILE_BASENAME}.${NUM_PARTICLES}.${NUM_RUNS}.dat)
