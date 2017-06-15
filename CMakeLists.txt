############################     BASE     ######################################

cmake_minimum_required(VERSION 2.8.11)
project(Nalu)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

############################     YAML     ######################################

find_library(YAML_LIBRARY NAMES yaml-cpp PATHS ${YAML_DIR}/lib )
find_path(YAML_INCLUDES yaml.h PATHS ${YAML_DIR}/include/yaml-cpp )
include_directories(${YAML_INCLUDES}/..)

############################     TRILINOS     ##################################

message("   Trilinos_DIR = ${Trilinos_DIR}")
set(CMAKE_PREFIX_PATH ${Trilinos_DIR} ${CMAKE_PREFIX_PATH})
set(CMAKE_PREFIX_PATH ${Trilinos_PATH} ${CMAKE_PREFIX_PATH})
find_package(Trilinos REQUIRED)

include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Trilinos_TPL_INCLUDE_DIRS})

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
set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})

# Enable Fortran if it is enabled in Trilinos
if(CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
endif()

# Build Nalu as shared libraries if Trilinos was compiled that way
if(Trilinos_BUILD_SHARED_LIBS)
  set(Nalu_LIBRARY_TYPE SHARED)
  set(BUILD_SHARED_LIBS ON)
  message("-- Building Nalu with shared libraries")
else(Trilinos_BUILD_SHARED_LIBS)
  set(Nalu_LIBRARY_TYPE STATIC)
endif(Trilinos_BUILD_SHARED_LIBS)

############################     NALU     ##################################

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(EXTRA_FLAGS "")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(EXTRA_FLAGS "")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(EXTRA_FLAGS "-restrict")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Trilinos_CXX_COMPILER_FLAGS} ${EXTRA_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Trilinos_Fortran_COMPILER_FLAGS}")
message("-- CMAKE_CXX_FLAGS     = ${CMAKE_CXX_FLAGS}")
message("-- CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")

file(GLOB SOURCE src/*.C src/*/*.C src/*/*.F)
file(GLOB HEADER include/*.h include/*/*.h)

include_directories(${CMAKE_SOURCE_DIR}/include)
add_library(nalu ${SOURCE} ${HEADER})
target_link_libraries(nalu ${Trilinos_LIBRARIES})
target_link_libraries(nalu ${YAML_LIBRARY})

file(GLOB UNIT_TESTS_SOURCES unit_tests/*.C unit_tests/*/*.C)

set(nalu_ex_name "naluX")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
   set(nalu_ex_name "naluXd")
   add_definitions("-Wall")
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

install(TARGETS ${utest_ex_name} ${nalu_ex_name} nalu
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)

############################     OPTIONS     ##################################

option(BUILD_DOCUMENTATION "Build documentation." OFF)
if(BUILD_DOCUMENTATION)
   add_subdirectory(docs)
endif()

option(ENABLE_TESTS "Enable regression testing." OFF)
if(ENABLE_TESTS)
   enable_testing()
   include(CTest)
   add_subdirectory(reg_tests)
endif()

MESSAGE("\nAnd CMake says...:")
