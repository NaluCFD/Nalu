# -----------------------------------------------------------
# -- Directories
# -----------------------------------------------------------

if(NOT ${NIGHTLY_DIR} STREQUAL "")
  message("Nightly test directory is ${NIGHTLY_DIR}")
else(NOT ${NIGHTLY_DIR} STREQUAL "")
  message( FATAL_ERROR "You need to set the NIGHTLY_DIR variable. CMake will exit." )
endif(NOT ${NIGHTLY_DIR} STREQUAL "")

set(NALU_DIR                        "${NIGHTLY_DIR}/Nalu")
set(NALURTEST_DIR                   "${NIGHTLY_DIR}/NaluRtest")
set(NALURTEST_MESH_DIR              "${NALURTEST_DIR}/mesh")
set(NALURTEST_XML_DIR               "${NALURTEST_DIR}/xml")
set(NALURTEST_NIGHTLY_DIR           "${NALURTEST_DIR}/nightly")
set(NALURTEST_PERF_DIR              "${NALURTEST_DIR}/performance")
set(TEST_RESULT_DIRECTORY           "${NIGHTLY_DIR}/runNaluRtest")
set(NIGHTLY_TEST_RESULT_DIRECTORY   "${TEST_RESULT_DIRECTORY}/nightly")
set(PERF_TEST_RESULT_DIRECTORY      "${TEST_RESULT_DIRECTORY}/performance")

# -----------------------------------------------------------
# -- REPOS
# -----------------------------------------------------------
set(NALU_REPO_URL                   "https://github.com/NaluCFD/Nalu.git")
set(NALURTEST_REPO_URL              "https://github.com/NaluCFD/NaluRtest.git")

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
set(CTEST_NALU_BINARY_NAME          "${CTEST_BINARY_DIRECTORY}/naluX")
set(CTEST_UNITTEST_BINARY_NAME      "${CTEST_BINARY_DIRECTORY}/unittestX")

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
if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone ${NALU_REPO_URL} ${CTEST_SOURCE_DIRECTORY}")
endif(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
       
## -- Update Command
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

## -- Configure Command
set(CTEST_CONFIGURE_COMMAND "cmake -DTrilinos_DIR:PATH=${TRILINOS_DIR} -DYAML_DIR:PATH=${YAML_DIR} -DENABLE_INSTALL:BOOL=OFF -DCMAKE_BUILD_TYPE=RELEASE ${CTEST_SOURCE_DIRECTORY}")

## -- Build Command
set(CTEST_BUILD_COMMAND "${MAKE} ${OPTION_BUILD}")

# -----------------------------------------------------------
# -- Configure CTest
# -----------------------------------------------------------

## -- Set TOLERANCE for test
set(TOLERANCE 0.000000001) # Default
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  if(${COMPILER_NAME} MATCHES "gcc")
    set(TOLERANCE 0.000000001)
  elseif(${COMPILER_NAME} MATCHES "clang")
    set(TOLERANCE 0.000000001)
  endif()
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  if(${COMPILER_NAME} MATCHES "gcc")
    set(TOLERANCE 0.000000001)
  elseif(${COMPILER_NAME} MATCHES "intel")
    set(TOLERANCE 0.000000001)
  endif()
endif()
message(" -- Using tolerance of ${TOLERANCE} --")

## -- CTest Test List and Customizations
configure_file(${CTEST_SOURCE_DIRECTORY}/cmake/CTestTestfile.cmake.in
               ${CTEST_BINARY_DIRECTORY}/CTestTestfile.cmake @ONLY)
file(COPY ${CTEST_SOURCE_DIRECTORY}/cmake/pass_fail.sh
     DESTINATION ${CTEST_BINARY_DIRECTORY}/)

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

## -- Clone (and pull) the test repo if necessary
if(NOT EXISTS "${NALURTEST_DIR}")
  execute_process(COMMAND "${CTEST_GIT_COMMAND}" "clone" "${NALURTEST_REPO_URL}" "${NALURTEST_DIR}"
                  WORKING_DIRECTORY ${NIGHTLY_DIR} )
endif(NOT EXISTS "${NALURTEST_DIR}")
execute_process(COMMAND "${CTEST_GIT_COMMAND}" "pull" WORKING_DIRECTORY ${NALURTEST_DIR})

## -- Prep test directory
message(" -- Prep test directory ${MODEL} - ${CTEST_BUILD_NAME} --")
include(${CTEST_SOURCE_DIRECTORY}/cmake/ctest_prepare_tests.cmake)

## -- Run CTest 
message(" -- Test ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_test(BUILD  "${CTEST_BINARY_DIRECTORY}"
           PARALLEL_LEVEL ${CTEST_PARALLEL_LEVEL}
           RETURN_VALUE res)

## -- Submit results to CDash
message(" -- Submit ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_submit( RETURN_VALUE res)

message(" -- Finished ${MODEL}  - ${CTEST_BUILD_NAME} --")
