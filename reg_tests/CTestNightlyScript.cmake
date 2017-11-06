if(NOT "${NIGHTLY_DIR}" STREQUAL "")
  message("Nightly test directory is ${NIGHTLY_DIR}")
else()
  message(FATAL_ERROR "You need to set the NIGHTLY_DIR variable. CMake will exit." )
endif()

if(NOT "${HOST_NAME}" STREQUAL "")
  message("Hostname is ${HOST_NAME}")
else()
  message(FATAL_ERROR "You need to set the HOST_NAME variable. CMake will exit." )
endif()

# -----------------------------------------------------------
# -- Configure CTest
# -----------------------------------------------------------

# Set important configuration variables
set(NALU_DIR "${NIGHTLY_DIR}/Nalu")
set(CTEST_SITE "${HOST_NAME}")
set(CTEST_BUILD_NAME "${CMAKE_SYSTEM_NAME}${EXTRA_BUILD_NAME}")
set(CTEST_SOURCE_DIRECTORY "${NALU_DIR}")
set(CTEST_BINARY_DIRECTORY "${NALU_DIR}/build")
set(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)
find_program(CTEST_GIT_COMMAND NAMES git)
find_program(MAKE NAMES make)

if("${BUILD_TYPE}" STREQUAL "")
  set(BUILD_TYPE "Release")
endif()

# Add parallelism capability to testing
include(ProcessorCount)
ProcessorCount(NP)
message(STATUS "\nNumber of processors detected: ${NP}")
set(CTEST_BUILD_FLAGS "-j${NP}")
set(CTEST_PARALLEL_LEVEL ${NP})

# Update Command
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")
set(CTEST_GIT_INIT_SUBMODULES "ON")

# Configure Command
set(CTEST_CONFIGURE_COMMAND "cmake -DTrilinos_DIR:PATH=${TRILINOS_DIR} -DYAML_DIR:PATH=${YAML_DIR} -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DENABLE_TESTS=ON ${TPL_TEST_ARGS} ${CTEST_SOURCE_DIRECTORY}")

# Build Command
set(CTEST_BUILD_COMMAND "${MAKE} ${CTEST_BUILD_FLAGS}")

# -----------------------------------------------------------
# -- Run CTest
# -----------------------------------------------------------

#ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
 
message("\n -- Start dashboard - ${CTEST_BUILD_NAME} --")
ctest_start("Nightly" TRACK "Nightly")

message("\n -- Update - ${CTEST_BUILD_NAME} --")
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE result)
message(" -- Update exit code = ${result} --")
if(result GREATER -1)
  message("\n -- Configure - ${CTEST_BUILD_NAME} --")
  ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE result)
  message(" -- Configure exit code = ${result} --")
  if(result EQUAL 0)
    message("\n -- Build - ${CTEST_BUILD_NAME} --")
    ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE result)
    message(" -- Build exit code = ${result} --")
    if(result EQUAL 0)
      # Need to have TMPDIR set to disk for building so it doesn't run out of space
      # but unset when running on these machines to stop OpenMPI from complaining
      string(COMPARE EQUAL "${HOST_NAME}" "peregrine.hpc.nrel.gov" is_equal_peregrine)
      string(COMPARE EQUAL "${HOST_NAME}" "merlin.hpc.nrel.gov" is_equal_merlin)
      if(is_equal_peregrine OR is_equal_merlin)
        message("Clearing TMPDIR variable...")
        unset(ENV{TMPDIR})
      endif()
      message("\n -- Test - ${CTEST_BUILD_NAME} --")
      ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}"
                 PARALLEL_LEVEL ${CTEST_PARALLEL_LEVEL}
                 RETURN_VALUE result)
      message(" -- Test exit code = ${result} --")
    endif()
  endif()
endif()

message("\n -- Submit - ${CTEST_BUILD_NAME} --")
set(CTEST_NOTES_FILES ${NIGHTLY_DIR}/jobs/nalu-test-log.txt)
ctest_submit(RETRY_COUNT 20
             RETRY_DELAY 20
             RETURN_VALUE result)
message(" -- Submit exit code = ${result} --")

message("\n -- Finished - ${CTEST_BUILD_NAME} --")
