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

SET(CMAKE_PREFIX_PATH ${TRILINOS_DIR} ${CMAKE_PREFIX_PATH})
FIND_PACKAGE(TRILINOS REQUIRED)

# This needs to be cleaned up
find_library(YAML_LIBRARY NAMES libyaml-cpp.a PATHS ${YAML_DIR}/lib )
find_path(YAML_INCLUDES yaml.h PATHS ${YAML_DIR}/include/yaml-cpp )
include_directories(${YAML_INCLUDES})

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

# ADD_OPTION(

# Optional Installation helpers
SET (INSTALL_NALU FALSE)
IF (ENABLE_INSTALL)
    SET (INSTALL_NALU TRUE)
ENDIF()

IF (INSTALL_NALU)
  set(BINARY_INSTALL_DIR bin)
  set(INCLUDE_INSTALL_DIR include)
  set(LIB_INSTALL_DIR lib)
  include(CMakePackageConfigHelpers)
ENDIF ()

MESSAGE("Setting and checking of compilers:")
SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

IF (CMAKE_Fortran_COMPILER) # Enable Fortran if it is enabled in Trilinos.
  ENABLE_LANGUAGE(Fortran)
ENDIF()

# Build Nalu as shared libraries if Trilinos was compiled that way

add_definitions(${Trilinos_CXX_COMPILER_FLAGS})
IF(Trilinos_BUILD_SHARED_LIBS)
  SET(Nalu_LIBRARY_TYPE SHARED)
  SET(BUILD_SHARED_LIBS ON)
  MESSAGE("-- Building Nalu with shared libraries")
ELSE(Trilinos_BUILD_SHARED_LIBS)
  SET(Nalu_LIBRARY_TYPE STATIC)
ENDIF(Trilinos_BUILD_SHARED_LIBS)

SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${Trilinos_CXX_COMPILER_FLAGS} )
MESSAGE("-- CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

#add_definitions(${Trilinos_CXX_COMPILER_FLAGS})

#######################     NALU     #######################

PROJECT(Nalu)

file (GLOB SOURCE src/*.C src/*/*.C src/*/*.F)
file (GLOB HEADER include/*.h include/*/*.h)

include_directories (${CMAKE_SOURCE_DIR}/include)
add_library (nalu ${SOURCE} ${HEADER})
target_link_libraries(nalu ${Trilinos_LIBRARIES})
target_link_libraries(nalu ${YAML_LIBRARY})

add_executable(naluX nalu.C)
target_link_libraries(naluX nalu)
MESSAGE("\nAnd CMake says...:")
