#=============================================================================
# Nightly tests
#=============================================================================

if(NOT ${NALURTEST_DIR} STREQUAL "")
  message("NaluRtest directory is ${NALURTEST_DIR}")
else(NOT ${NALURTEST_DIR} STREQUAL "")
  message( FATAL_ERROR "You need to set the ${NALURTEST_DIR} variable. CMake will exit." )
endif()

if(NOT ${RUNNALURTEST_DIR} STREQUAL "")
  message("runNaluRtest directory is ${RUNNALURTEST_DIR}")
else(NOT ${RUNNALURTEST_DIR} STREQUAL "")
  message( FATAL_ERROR "You need to set the ${RUNNALURTEST_DIR} variable. CMake will exit." )
endif()

set(NALURTEST_MESH_DIR              "${NALURTEST_DIR}/mesh")
set(NALURTEST_XML_DIR               "${NALURTEST_DIR}/xml")
set(NALURTEST_NIGHTLY_DIR           "${NALURTEST_DIR}/nightly")
set(NALURTEST_PERF_DIR              "${NALURTEST_DIR}/performance")
set(TEST_RESULT_DIRECTORY           "${RUNNALURTEST_DIR}")
set(NIGHTLY_TEST_RESULT_DIRECTORY   "${TEST_RESULT_DIRECTORY}/nightly")
set(PERF_TEST_RESULT_DIRECTORY      "${TEST_RESULT_DIRECTORY}/performance")

## -- Set TOLERANCE for test
set(TOLERANCE 0.0000001) # Default
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  if(${COMPILER_NAME} MATCHES "gcc")
    set(TOLERANCE 0.0000001)
  elseif(${COMPILER_NAME} MATCHES "clang")
    set(TOLERANCE 0.0000001)
  endif()
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  if(${COMPILER_NAME} MATCHES "gcc")
    set(TOLERANCE 0.0000001)
  elseif(${COMPILER_NAME} MATCHES "intel")
    set(TOLERANCE 0.0001)
  endif()
endif()
message(" -- Using tolerance of ${TOLERANCE} --")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake.in
               ${CMAKE_BINARY_DIR}/CTestCustom.cmake @ONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestPrepareTests.cmake.in
               ${CMAKE_BINARY_DIR}/CTestPrepareTests.cmake @ONLY)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh
     DESTINATION ${CMAKE_BINARY_DIR}/)

include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
