if(NOT ${NIGHTLY_DIR} STREQUAL "")
  message("Nightly test directory is ${NIGHTLY_DIR}")
else()
  message(FATAL_ERROR "You need to set the NIGHTLY_DIR variable. CMake will exit." )
endif()

if(NOT ${HOST_NAME} STREQUAL "")
  message("Hostname is ${HOST_NAME}")
else()
  message(FATAL_ERROR "You need to set the HOST_NAME variable. CMake will exit." )
endif()

# -----------------------------------------------------------
# -- Configure CTest
# -----------------------------------------------------------

# Set important configuration variables
set(MODEL "nightly")
set(NALU_DIR  "${NIGHTLY_DIR}/Nalu")
set(CTEST_SITE "${HOST_NAME}")
set(CTEST_BUILD_NAME "${CMAKE_SYSTEM_NAME}${EXTRA_BUILD_NAME}")
set(CTEST_SOURCE_DIRECTORY "${NALU_DIR}")
set(CTEST_BINARY_DIRECTORY "${NALU_DIR}/build")
find_program(CTEST_GIT_COMMAND NAMES git)
find_program(MAKE NAMES make)

# Add parallelism capability to testing
include(ProcessorCount)
ProcessorCount(NP)
message(STATUS "Number of processors detected: ${NP}")
set(CTEST_BUILD_FLAGS -j${NP})
set(CTEST_PARALLEL_LEVEL ${NP})
set(MAKE_OPTIONS "-j${NP}")

# Update Command
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")
set(CTEST_GIT_INIT_SUBMODULES "ON")

# Configure Command
set(CTEST_CONFIGURE_COMMAND "cmake -DTrilinos_DIR:PATH=${TRILINOS_DIR} -DYAML_DIR:PATH=${YAML_DIR} -DENABLE_INSTALL:BOOL=OFF -DCMAKE_BUILD_TYPE=RELEASE -DENABLE_TESTS=ON ${CTEST_SOURCE_DIRECTORY}")

# Build Command
set(CTEST_BUILD_COMMAND "${MAKE} ${MAKE_OPTIONS}")

# -----------------------------------------------------------
# -- Run CTest
# -----------------------------------------------------------

ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
 
message(" -- Start dashboard ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_start(${MODEL} TRACK ${MODEL})

message(" -- Update ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)

message(" -- Configure ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)

message(" -- Build ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)

message(" -- Test ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}"
           PARALLEL_LEVEL ${CTEST_PARALLEL_LEVEL}
           RETURN_VALUE res)

message(" -- Submit ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_submit(RETRY_COUNT 20
             RETRY_DELAY 20
             RETURN_VALUE res)

message(" -- Finished ${MODEL} - ${CTEST_BUILD_NAME} --")
