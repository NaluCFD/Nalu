
#	This is the top level CMakeLists.txt file 
#	It should be located within master/nalu
#	But not in the BUILD folder itself
#
###########################################################################################################

cmake_minimum_required (VERSION 2.8.8)

###########################################################################################################
#
#  THESE PATHS MAY NEED TO BE CHANGED DEPENDING ON WHERE THEY ARE STORED!!! 
#

SET(CMAKE_PREFIX_PATH "/usr/local/Cellar/cmake" ${TRILINOS_DIR} ${CMAKE_PREFIX_PATH}) 
FIND_PACKAGE(TRILINOS PATHS </Users/nalu/install/lib/cmake>)

SET(CMAKE_FIND_LIBRARY_PREFIXES lib)
SET(CMAKE_FIND_LIBRARY_SUFFIXES a) 

find_library(YAML_LIBRARY NAMES libyaml-cpp.a PATHS /Users/nalu/install/lib )
find_path(YAML_INCLUDES yaml.h PATHS /Users/nalu/install/include/yaml-cpp )  
include_directories(${YAML_INCLUDES}) 

include(FindPackageHandleStandardArgs)

include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(/Users/nalu/install/include/yaml-cpp)
include_directories(/usr/local/Cellar/superlu/4.3/include/superlu) 
include_directories(${CMAKE_SOURCE_DIR}/../utility/include) 
include_directories(${CMAKE_SOURCE_DIR}/../TPLs_src/Trilinos/packages/stk/stk_transfer)

#######################################     TRILINOS     ##################################################

MESSAGE("\nFound Trilinos!  Here are the details: ")
# Uncomment lines to view more details about Trilinos 
#MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
#MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
#MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
#MESSAGE("   Trilinos_BIN_DIRS = ${Trilinos_BIN_DIRS}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
#MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
#MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
MESSAGE("End of Trilinos details\n")


# Optional Installation helpers
SET (INSTALL_NALU FALSE)
IF (ENABLE_INSTALL)
  IF (${CMAKE_VERSION} VERSION_GREATER 2.8.8 OR ${CMAKE_VERSION} VERSION_EQUAL 2.8.8)
    SET (INSTALL_NALU TRUE)
    MESSAGE("")
    MESSAGE("Installation files for Nalu will be created")
    MESSAGE("-- install location = ${CMAKE_INSTALL_PREFIX}")
    MESSAGE("")
  ELSE()
    MESSAGE("")
    MESSAGE("Installation has been requested, but cannot continue")
    MESSAGE("-- Please upgrade CMAKE to version >= 2.8.8")
    MESSAGE("")
  ENDIF()
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

IF(Trilinos_BUILD_SHARED_LIBS)
  SET(Nalu_LIBRARY_TYPE SHARED)
  SET(BUILD_SHARED_LIBS ON)
  MESSAGE("-- Building Nalu with shared libraries")
ELSE(Trilinos_BUILD_SHARED_LIBS)
  SET(Nalu_LIBRARY_TYPE STATIC)
ENDIF(Trilinos_BUILD_SHARED_LIBS)

IF (NOT CMAKE_CXX_FLAGS)
  SET(CMAKE_CXX_FLAGS ${Trilinos_CXX_COMPILER_FLAGS} )
ENDIF()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NALU_CXX_FLAGS}")
MESSAGE("   CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
add_definitions(${CMAKE_CXX_FLAGS})


#######################################     NALU     ######################################################

PROJECT(Nalu)

file (GLOB SOURCE src/*.C src/*/*.C src/*/*.F)
file (GLOB HEADER include/*.h include/*/*.h)

include_directories (${CMAKE_SOURCE_DIR}/include)
add_library (nalu ${SOURCE} ${HEADER})
target_link_libraries(nalu ${Trilinos_LIBRARIES})

add_executable(nalu_exe nalu.C)
target_link_libraries(nalu_exe nalu)
target_link_libraries(nalu_exe ${YAML_LIBRARY}) 
#target_link_libraries(nalu_exe ${STK_MESH_LIBRARY})
MESSAGE("\nAnd CMake says...:")
