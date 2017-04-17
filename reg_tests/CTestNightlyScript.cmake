# -----------------------------------------------------------
# -- Directories
# -----------------------------------------------------------

if(NOT ${NIGHTLY_DIR} STREQUAL "")
  message("Nightly test directory is ${NIGHTLY_DIR}")
else(NOT ${NIGHTLY_DIR} STREQUAL "")
  message( FATAL_ERROR "You need to set the NIGHTLY_DIR variable. CMake will exit." )
endif(NOT ${NIGHTLY_DIR} STREQUAL "")

set(NALU_DIR                        "${NIGHTLY_DIR}/Nalu")

# -----------------------------------------------------------
# -- REPOS
# -----------------------------------------------------------
set(NALU_REPO_URL                   "https://github.com/NaluCFD/Nalu.git")

# -----------------------------------------------------------
# -- Get environment
# -----------------------------------------------------------

## -- Set hostname
find_program(HOSTNAME_CMD NAMES hostname)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  exec_program("${HOSTNAME_CMD} -f" ARGS OUTPUT_VARIABLE HOSTNAME)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  exec_program("${HOSTNAME_CMD} -d" ARGS OUTPUT_VARIABLE HOSTNAME)
endif()
if(${HOSTNAME} MATCHES "hpc.nrel.gov")
  set(HOSTNAME "peregrine.hpc.nrel.gov")
endif()

## -- Set site / build name
find_program(UNAME NAMES uname)
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(getuname)
getuname(CPU -m)

set(CTEST_SITE "${HOSTNAME}")
set(CTEST_BUILD_NAME "${CMAKE_SYSTEM_NAME}-${CPU}-${COMPILER_NAME}-${TRILINOS_BRANCH}")

## -- Git command
find_program(CTEST_GIT_COMMAND NAMES git)

## -- make command
find_program(MAKE NAMES make)

# -----------------------------------------------------------
# -- build specific
# -----------------------------------------------------------

set(MODEL                           "nightly")
set(CTEST_SOURCE_DIRECTORY          "${NALU_DIR}")
set(CTEST_BINARY_DIRECTORY          "${NALU_DIR}/build")

## -- Build options 
include(ProcessorCount)
ProcessorCount(NP)
message(STATUS "Number of processors detected: ${NP}")
set(CTEST_BUILD_FLAGS               -j${NP})
set(CTEST_PARALLEL_LEVEL            ${NP})
set(OPTION_BUILD                    "-j${NP}")

# -----------------------------------------------------------
# -- commands
# -----------------------------------------------------------

## -- Checkout command
#if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
#  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone ${NALU_REPO_URL} ${CTEST_SOURCE_DIRECTORY}")
#endif(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
       
## -- Update Command
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")
set(CTEST_GIT_INIT_SUBMODULES "ON")

## -- Configure Command
set(CTEST_CONFIGURE_COMMAND "cmake -DTrilinos_DIR:PATH=${TRILINOS_DIR} -DYAML_DIR:PATH=${YAML_DIR} -DENABLE_INSTALL:BOOL=OFF -DCMAKE_BUILD_TYPE=RELEASE -DENABLE_TESTS=ON ${CTEST_SOURCE_DIRECTORY}")

## -- Build Command
set(CTEST_BUILD_COMMAND "${MAKE} ${OPTION_BUILD}")

# -----------------------------------------------------------
# -- Run CTest
# -----------------------------------------------------------

## -- Start
message(" -- Start dashboard ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_start(${MODEL} TRACK ${MODEL})

## -- Update
message(" -- Update ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)

## -- Configure
message(" -- Configure ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_configure(BUILD  "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)

## -- Build
message(" -- Build ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_build(BUILD  "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)

## -- Run CTest 
message(" -- Test ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_test(BUILD  "${CTEST_BINARY_DIRECTORY}"
           PARALLEL_LEVEL ${CTEST_PARALLEL_LEVEL}
           RETURN_VALUE res)

## -- Submit results to CDash
message(" -- Submit ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_submit(RETRY_COUNT 20
             RETRY_DELAY 20
             RETURN_VALUE res)

message(" -- Finished ${MODEL}  - ${CTEST_BUILD_NAME} --")
