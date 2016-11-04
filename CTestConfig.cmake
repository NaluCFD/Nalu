
set(CTEST_PROJECT_NAME "Nalu")

# Time at which to run the test
set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")

# Set these according to your CDash server
# make sure you have the following in CMakeLists.txt
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_DROP_METHOD http)
set(CTEST_DROP_SITE "my.cdash.org")
set(CTEST_DROP_LOCATION "/submit.php?project=Nalu")
set(CTEST_DROP_SITE_CDASH TRUE)

# Get a parallel build process
include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
 set(CTEST_BUILD_FLAGS -j${N})
 set(CTEST_PARALLEL_LEVEL ${N})
 message(STATUS "Number of Processors detected: ${N}")
endif()