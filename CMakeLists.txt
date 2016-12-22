#
#	This is the top level CMakeLists.txt file
#	It should be located within master/nalu
#	But not in the BUILD folder itself
#
###########################################################################################################

cmake_minimum_required (VERSION 2.8.11)

###########################################################################################################
#
#  THESE PATHS MAY NEED TO BE CHANGED DEPENDING ON WHERE THEY ARE STORED!!!
#

SET(CMAKE_FIND_LIBRARY_PREFIXES lib)
SET(CMAKE_FIND_LIBRARY_SUFFIXES a)

MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
SET(CMAKE_PREFIX_PATH ${Trilinos_DIR} ${CMAKE_PREFIX_PATH})
SET(CMAKE_PREFIX_PATH ${Trilinos_PATH} ${CMAKE_PREFIX_PATH})
find_package(Trilinos REQUIRED)

# This needs to be cleaned up
find_library(YAML_LIBRARY NAMES libyaml-cpp.a PATHS ${YAML_DIR}/lib )
find_path(YAML_INCLUDES yaml.h PATHS ${YAML_DIR}/include/yaml-cpp )
include_directories(${YAML_INCLUDES}/..)

#include(FindPackageHandleStandardArgs)

include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Trilinos_TPL_INCLUDE_DIRS})

#######################################     TRILINOS	 ##################################################

MESSAGE("\nFound Trilinos!  Here are the details: ")
# Uncomment lines to view more details about Trilinos
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_BIN_DIRS = ${Trilinos_BIN_DIRS}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
MESSAGE("End of Trilinos details\n")

# Optional Installation helpers
SET (INSTALL_NALU FALSE)
IF (ENABLE_INSTALL)
    SET (INSTALL_NALU TRUE)
ENDIF()

MESSAGE("Setting and checking of compilers:")
SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

IF (CMAKE_Fortran_COMPILER) # Enable Fortran if it is enabled in Trilinos.
  ENABLE_LANGUAGE(Fortran)
ENDIF()

# Build Nalu as shared libraries if Trilinos was compiled that way

IF(Trilinos_BUILD_SHARED_LIBS)
  SET(Nalu_LIBRARY_TYPE SHARED)
  SET(BUILD_SHARED_LIBS ON)
  MESSAGE("-- Building Nalu with shared libraries")
ELSE(Trilinos_BUILD_SHARED_LIBS)
  SET(Nalu_LIBRARY_TYPE STATIC)
ENDIF(Trilinos_BUILD_SHARED_LIBS)

#######################################     NALU     ##########################

PROJECT(Nalu)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Trilinos_CXX_COMPILER_FLAGS}")
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Trilinos_Fortran_COMPILER_FLAGS}")
MESSAGE("-- CMAKE_CXX_FLAGS     = ${CMAKE_CXX_FLAGS}")
MESSAGE("-- CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")

file (GLOB SOURCE src/*.C src/*/*.C src/*/*.F)
file (GLOB HEADER include/*.h include/*/*.h)

include_directories (${CMAKE_SOURCE_DIR}/include)
add_library (nalu ${SOURCE} ${HEADER})
target_link_libraries(nalu ${Trilinos_LIBRARIES})
target_link_libraries(nalu ${YAML_LIBRARY})

file (GLOB UNIT_TESTS_SOURCES unit_tests/*.C)

set(nalu_ex_name "naluX")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
   set (nalu_ex_name "naluXd")
   add_definitions("-Wall")
   message("Debug Build")
endif()

add_executable(${nalu_ex_name} nalu.C)
target_link_libraries(${nalu_ex_name} nalu)

set(utest_ex_name "unittestX")
if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
   set(utest_ex_name "unittestXd")
endif()

add_executable(${utest_ex_name} unit_tests.C ${UNIT_TESTS_SOURCES})
target_link_libraries(${utest_ex_name} nalu)

IF (INSTALL_NALU)
  install(TARGETS ${utest_ex_name} ${nalu_ex_name} nalu
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)
  install(DIRECTORY include/ DESTINATION include)
  include(CMakePackageConfigHelpers)
ENDIF ()

MESSAGE("\nAnd CMake says...:")
