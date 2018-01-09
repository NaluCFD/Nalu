if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "-w ${CMAKE_Fortran_FLAGS}")
endif()
add_library(DISCON SHARED DISCON.f90)
