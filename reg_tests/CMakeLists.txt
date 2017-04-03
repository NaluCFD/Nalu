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

## -- Set TOLERANCE for testing
if(NOT ${TEST_TOLERANCE} STREQUAL "")
  set(TOLERANCE ${TEST_TOLERANCE}) # User defined
else(NOT ${TEST_TOLERANCE} STREQUAL "")
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
      set(TOLERANCE 0.0000001)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      set(TOLERANCE 0.0000001)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      set(TOLERANCE 0.0000001)
    else()
      set(TOLERANCE 0.0000001) # Mac default
    endif()
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
      set(TOLERANCE 0.0000001)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      set(TOLERANCE 0.0000000000001)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      set(TOLERANCE 0.0001)
    else()
      set(TOLERANCE 0.0000001) # Linux default
    endif()
  endif()
endif()
message(" -- Using test tolerance of ${TOLERANCE} --")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake.in
               ${CMAKE_BINARY_DIR}/CTestCustom.cmake @ONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestPrepareTests.cmake.in
               ${CMAKE_BINARY_DIR}/CTestPrepareTests.cmake @ONLY)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh
     DESTINATION ${CMAKE_BINARY_DIR}/)

include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
