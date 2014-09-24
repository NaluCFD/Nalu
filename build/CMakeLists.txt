cmake_minimum_required (VERSION 2.8.8)

# Get Trilinos as one entity
SET(CMAKE_PREFIX_PATH ${TRILINOS_DIR} ${CMAKE_PREFIX_PATH})
FIND_PACKAGE(Trilinos REQUIRED PATHS /usr/local/Cellar/Trilinos)



MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
#MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
#MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_BIN_DIRS = ${Trilinos_BIN_DIRS}")
#MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
#MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
#MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
#MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
MESSAGE("End of Trilinos details\n")


# Optional Installation helpers
# Note that some of this functionality depends on CMAKE > 2.8.8
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



###########################################################################################################
PROJECT(Nalu)


IF(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
  SET(USING_INTEL TRUE)
# On Intel, turn of interprocedural optimization at linking to keep linker from hanging
# This could possibly be removed with further study (GAH)
#  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -no-ip")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd2536,2196,858,47")
ENDIF()

MESSAGE("\nChecking Trilinos build for required and optional packages:")


### Check for some OPTIONAL trilinos packages, which set DEFINES in Alany build ###

# MPI check -- defaults to TRUE
# MESSAGE("-- Checking if MPI is enabled in Trilinos:")
LIST(FIND Trilinos_TPL_LIST MPI MPI_List_ID)
IF (MPI_List_ID GREATER -1)
  MESSAGE("-- Looking for MPI in Trilinos:    Found, compiling with -DALBANY_MPI")
  SET(NALU_MPI TRUE)
  ADD_DEFINITIONS(-DNALU_MPI)
ELSE()
  MESSAGE("-- Looking for MPI in Trilinos:    NOT found, MPI NOT ENABLED!")
  SET(NALU_MPI FALSE)
ENDIF()

# set optional dependency on stk, for uniform adaptivity
LIST(FIND Trilinos_PACKAGE_LIST STK stk_List_ID)
IF (stk_List_ID GREATER -1)
  ADD_DEFINITIONS(-DNALU_STK)
  MESSAGE("-- Looking for stk:        Found, compiling with -DNALU_STK")
  SET(NALU_STK TRUE)
ELSE()
  MESSAGE("-- Looking for stk:        NOT found.")
  SET(NALU_STK FALSE)
ENDIF()

# set optional dependency on SEACAS, for I/O via stk_io
# MESSAGE("-- Looking for SEACAS:")
#LIST(FIND Trilinos_PACKAGE_LIST SEACASIoss SEACASIoss_List_ID)
#IF (SEACASIoss_List_ID GREATER -1)
#  ADD_DEFINITIONS(-DNALU_SEACAS)
#  MESSAGE("-- Looking for SEACASIoss:         Found, compiling with -DALBANY_SEACAS")
#  SET(NALU_SEACAS TRUE)
#ELSE()
#  MESSAGE("-- Looking for SEACASIoss:         NOT found.")
#  SET(NALU_SEACAS FALSE)
#ENDIF()


file (GLOB SOURCE src/*.C)
file (GLOB HEADER include/*.h)

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/src)
add_library (nalu ${SOURCE} ${HEADER})




# Make a clean build target for convenience

ADD_CUSTOM_TARGET( cmake_clean
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/Makefile
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfile.cmake
)

IF(NALU_MPI)
  IF(NALU_MPI_OPTIONS)
    SET(MPIEX  ${NALU_MPI_EXEC})
    SET(MPINPF ${NALU_MPI_EXEC_NUMPROCS_FLAG})
    SET(MPIMNP ${NALU_MPI_EXEC_MAX_NUMPROCS})
    SET(MPIPRE ${NALU_MPI_LEADING_OPTIONS})
    # For now - assume NALU_MPI_OPTIONS means to ONLY attempt the parallel
    # tests
    SET(NALU_PARALLEL_ONLY true)
    MESSAGE("-- CTest mpi_exec call:  ${MPIEX} ${MPIPRE} ${MPINPF} ${MPIMNP}")
  ELSE()
    SET(MPIEX  ${Trilinos_MPI_EXEC})
    SET(MPINPF ${Trilinos_MPI_EXEC_NUMPROCS_FLAG})
    SET(MPIMNP ${Trilinos_MPI_EXEC_MAX_NUMPROCS})
    SET(MPIPRE "")
    SET(NALU_PARALLEL_ONLY false)
    MESSAGE("-- CTest mpi_exec call:  ${MPIEX} ${MPINPF} ${MPIMNP}")
  ENDIF()
ENDIF(NALU_MPI)

MESSAGE("\nAnd CMake says...:")
