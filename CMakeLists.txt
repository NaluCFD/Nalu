############################ BASE ######################################
cmake_minimum_required (VERSION 3.23.0)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
cmake_policy(SET CMP0057 NEW)

########################## OPTIONS #####################################
option(ENABLE_TESTS "Enable regression testing." OFF)
option(ENABLE_DOCUMENTATION "Build documentation." OFF)
option(ENABLE_SPHINX_API_DOCS "Link Doxygen API docs to Sphinx" OFF)
option(ENABLE_PARAVIEW_CATALYST
      "Enable ParaView Catalyst. Requires external installation of Trilinos Catalyst IOSS adapter."
       OFF)
option(ENABLE_WARNINGS "Add -Wall to show compiler warnings" ON)
option(ENABLE_EXTRA_WARNINGS "Add -Wextra to show even more compiler warnings" OFF)

########################### NALU-first #####################################
# Set Nalu's compilers, CMAKE_FIND_LIBRARY_PREFIXES
# and CMAKE_FIND_LIBRARY_SUFFIXES, etc, after inheriting settings from Trilinos
project(Nalu CXX Fortran)

########################## TRILINOS ####################################
message("   Trilinos_DIR = ${Trilinos_DIR}")
set(CMAKE_PREFIX_PATH ${Trilinos_DIR} ${CMAKE_PREFIX_PATH})
set(CMAKE_PREFIX_PATH ${Trilinos_PATH} ${CMAKE_PREFIX_PATH})
find_package(Trilinos NO_MODULE NO_CMAKE_ENVIRONMENT_PATH NO_SYSTEM_ENVIRONMENT_PATH REQUIRED)
include_directories(SYSTEM ${Trilinos_INCLUDE_DIRS})
include_directories(SYSTEM ${Trilinos_TPL_INCLUDE_DIRS})

message("\nFound Trilinos!  Here are the details: ")
message("   Trilinos_DIR = ${Trilinos_DIR}")
message("   Trilinos_VERSION = ${Trilinos_VERSION}")
message("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
message("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
message("   Trilinos_BIN_DIRS = ${Trilinos_BIN_DIRS}")
message("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
message("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
message("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
message("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
message("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
message("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
message("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
message("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
message("End of Trilinos details\n")

message("Setting and checking of compilers:")
set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

# Build Nalu as shared libraries if that is how Trilinos was built
if(Trilinos_BUILD_SHARED_LIBS)
  set(Nalu_LIBRARY_TYPE SHARED)
  set(BUILD_SHARED_LIBS ON)
  message("-- Building Nalu with shared libraries")
else(Trilinos_BUILD_SHARED_LIBS)
  set(Nalu_LIBRARY_TYPE STATIC)
endif(Trilinos_BUILD_SHARED_LIBS)

############################ YAML ######################################
set(CMAKE_PREFIX_PATH ${YAML_DIR} ${CMAKE_PREFIX_PATH})
find_package(YAML-CPP QUIET)
if(YAML-CPP_FOUND)
  # YAML master branch is used
  include_directories(SYSTEM ${YAML_CPP_INCLUDE_DIR})
else()
  # YAML 0.5.3 is used
  find_library(YAML_CPP_LIBRARIES NAMES yaml-cpp PATHS ${YAML_DIR}/lib)
  find_path(YAML_CPP_INCLUDE_DIR yaml.h PATHS ${YAML_DIR}/include/yaml-cpp)
  if((DEFINED YAML_CPP_LIBRARIES) AND (DEFINED YAML_CPP_INCLUDE_DIR))
    include_directories(SYSTEM ${YAML_CPP_INCLUDE_DIR}/..)
    set(YAML-CPP_FOUND TRUE)
  endif()
endif()
if(YAML-CPP_FOUND)
  message("-- Found YAML-CPP = ${YAML_DIR}")
else()
  message(FATAL_ERROR "YAML-CPP NOT FOUND")
endif()

########################### NALU-second #####################################

# Add any extra flags based on compiler and/or OS
message("-- CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message("-- CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  set(EXTRA_CXX_FLAGS "")
  set(EXTRA_Fortran_FLAGS "")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(EXTRA_CXX_FLAGS "")
  set(EXTRA_Fortran_FLAGS "")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(EXTRA_CXX_FLAGS "")
  set(EXTRA_Fortran_FLAGS "")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(EXTRA_CXX_FLAGS "-restrict")
  set(EXTRA_Fortran_FLAGS "")
endif()

if(ENABLE_WARNINGS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall")
endif()
if(ENABLE_EXTRA_WARNINGS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -pedantic")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wextra -pedantic")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Trilinos_CXX_COMPILER_FLAGS} ${EXTRA_CXX_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Trilinos_Fortran_COMPILER_FLAGS} ${EXTRA_Fortran_FLAGS}")

message("-- CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message("-- CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")

################################# TAG VERSIONS #####################################
include(GetGitRevisionDescription)
get_git_head_revision(NALU_GITREFSPEC NALU_GIT_COMMIT_SHA)
if (NALU_GIT_COMMIT_SHA)
  git_describe(NALU_VERSION_TAG "--tags")
  git_local_changes(NALU_REPO_DIRTY)
  option(NALU_HAVE_GIT_INFO "Git version tagging for Nalu" ON)
  if (${NALU_VERSION_TAG} MATCHES ".*-NOTFOUND")
    set(NALU_VERSION_TAG "v1.2.0")
  endif()
endif()

find_file(TRILINOS_HAVE_GIT_INFO
  "TrilinosRepoVersion.txt"
  PATHS "${Trilinos_DIR}/../../../"
  DOC "Check if Trilinos Git version info exists"
  NO_DEFAULT_PATH)
if (NOT ${TRILINOS_HAVE_GIT_INFO} MATCHES ".*-NOTFOUND")
  file(STRINGS "${Trilinos_DIR}/../../../TrilinosRepoVersion.txt"
    TRILINOS_REPO_VERSION_TXT LIMIT_INPUT 1024)
  list(GET TRILINOS_REPO_VERSION_TXT 1 TRILINOS_REPO_COMMIT_STR)
  string(REGEX MATCH "^[a-z0-9]+" TRILINOS_GIT_COMMIT_SHA ${TRILINOS_REPO_COMMIT_STR})
  set(TRILINOS_VERSION_TAG "${Trilinos_VERSION}-g${TRILINOS_GIT_COMMIT_SHA}")
else()
  set(TRILINOS_VERSION_TAG "${Trilinos_VERSION}")
endif()
string(TIMESTAMP NALU_VERSION_TIMESTAMP "%Y-%m-%d %H:%M:%S (UTC)" UTC)
configure_file("${CMAKE_SOURCE_DIR}/cmake/NaluVersionInfo.h.in"
  "${CMAKE_BINARY_DIR}/include/NaluVersionInfo.h" @ONLY)
#### END TAG VERSIONS

file(GLOB SOURCE src/*.C src/*/*.C src/*/*.F)
file(GLOB HEADER include/*.h include/*/*.h
  ${CMAKE_CURRENT_BINARY_DIR}/include/*.h)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/include)
add_library(nalu ${SOURCE} ${HEADER})
target_link_libraries(nalu ${Trilinos_LIBRARIES})
target_link_libraries(nalu ${YAML_CPP_LIBRARIES})

file(GLOB UNIT_TESTS_SOURCES unit_tests/*.C unit_tests/*/*.C)

set(nalu_ex_name "naluX")
message("-- CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
   set(nalu_ex_name "naluXd")
   message("Debug Build")
endif()

add_executable(${nalu_ex_name} nalu.C)
target_link_libraries(${nalu_ex_name} nalu)

set(utest_ex_name "unittestX")
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
   set(utest_ex_name "unittestXd")
endif()

add_executable(${utest_ex_name} unit_tests.C ${UNIT_TESTS_SOURCES})
target_link_libraries(${utest_ex_name} nalu)
target_include_directories(${utest_ex_name} PUBLIC "${CMAKE_SOURCE_DIR}/unit_tests")

set(nalu_ex_catalyst_name "naluXCatalyst")
if(ENABLE_PARAVIEW_CATALYST)
   set(PARAVIEW_CATALYST_INSTALL_PATH
       ""
       CACHE
       PATH
       "Path to external installation of Trilinos Catalyst IOSS plugin.")
   configure_file(cmake/naluXCatalyst.in ${nalu_ex_catalyst_name} @ONLY)
  add_definitions("-DNALU_USES_CATALYST")
endif()

install(TARGETS ${utest_ex_name} ${nalu_ex_name} nalu
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)
if(ENABLE_PARAVIEW_CATALYST)
  install(PROGRAMS ${CMAKE_BINARY_DIR}/naluXCatalyst DESTINATION bin)
endif()

include(CMakePackageConfigHelpers)

if(ENABLE_DOCUMENTATION)
   add_subdirectory(docs)
endif()

if(ENABLE_TESTS)
   message(" Enabling regression test suite")
   enable_testing()
   include(CTest)
   add_subdirectory(reg_tests)
endif()

message("\nAnd CMake says...:")
