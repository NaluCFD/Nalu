# Find HYPRE linear solver library
#
# Set HYPRE_DIR to the base directory where the package is installed
#
# Sets two variables
#   - HYPRE_INCLUDE_DIRS
#   - HYPRE_LIBRARIES
#

find_path(HYPRE_INCLUDE_DIRS
  HYPRE.h
  HINTS ${HYPRE_DIR} ${CMAKE_INSTALL_PREFIX}
  PATH_SUFFIXES include)

find_library(HYPRE_LIBRARIES
  NAMES HYPRE HYPRE
  HINTS ${HYPRE_DIR} ${CMAKE_INSTALL_PREFIX}
  PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  HYPRE DEFAULT_MSG HYPRE_INCLUDE_DIRS HYPRE_LIBRARIES)
mark_as_advanced(HYPRE_INCLUDE_DIRS HYPRE_LIBRARIES)
