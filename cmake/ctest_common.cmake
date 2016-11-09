# -----------------------------------------------------------
# -- Directories
# -----------------------------------------------------------

set(NIGHTLY_DIR                      "/scratch/mhenryde/nightlyBuildAndTest")
set(NALU_DIR                         "${NIGHTLY_DIR}/Nalu")
set(NALURTEST_DIR                    "${NIGHTLY_DIR}/NaluRtest")
set(NIGHTLY_TEST_RESULT_DIRECTORY    "${NIGHTLY_DIR}/runNaluRtest/nightly")


# -----------------------------------------------------------
# -- VARIABLES
# -----------------------------------------------------------
set(NALU_REPO_URL                       "git@github.com:marchdf/Nalu.git")
set(NALURTEST_REPO_URL                  "https://github.com/NaluCFD/NaluRtest.git")

# -----------------------------------------------------------
# -- Get environment
# -----------------------------------------------------------

## -- Set hostname
## --------------------------
find_program(HOSTNAME_CMD NAMES hostname)
exec_program(${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME)

set(CTEST_SITE                          "${HOSTNAME}")

## -- Set site / build name
## --------------------------
find_program(UNAME NAMES uname)
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(getuname)
  
getuname(osname -s)
getuname(osrel  -r)
getuname(cpu    -m)

set(CTEST_BUILD_NAME                    "${osname}-${cpu}-${compiler}")

## -- Git command
## ----------------
find_program(CTEST_GIT_COMMAND NAMES git)

## -- make command
## -----------------
find_program(MAKE NAMES make)

# -----------------------------------------------------------
# -- build specific
# -----------------------------------------------------------

set(MODEL                               "nightly")

## -- SRC Dir
set(CTEST_SOURCE_DIRECTORY              "${NALU_DIR}")

## -- BIN Dir
set(CTEST_BINARY_DIRECTORY              "${NALU_DIR}/build")

## -- Build options 
include(ProcessorCount)
ProcessorCount(NP)
message(STATUS "Number of processors detected: ${NP}")
set(CTEST_BUILD_FLAGS                   -j${NP})
set(CTEST_PARALLEL_LEVEL                ${NP})
set(OPTION_BUILD                        "-j${NP}")

# -----------------------------------------------------------
# -- commands
# -----------------------------------------------------------

## -- Checkout command
if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
        set(CTEST_CHECKOUT_COMMAND     "${CTEST_GIT_COMMAND} clone ${NALU_REPO_URL} ${CTEST_SOURCE_DIRECTORY}")
endif(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
       
## -- Update Command
set(CTEST_UPDATE_COMMAND               "${CTEST_GIT_COMMAND}")

## -- Configure Command
set(CTEST_CONFIGURE_COMMAND            "${CTEST_BINARY_DIRECTORY}/do-configNalu_release-${compiler}")

## -- Build Command
set(CTEST_BUILD_COMMAND                "${MAKE} ${OPTION_BUILD}")

# -----------------------------------------------------------
# -- Configure CTest
# -----------------------------------------------------------

## -- CTest Testfile
configure_file(${CTEST_SOURCE_DIRECTORY}/cmake/CTestTestfile.cmake ${CTEST_BINARY_DIRECTORY}/CTestTestfile.cmake)

# -----------------------------------------------------------
# -- Start CTest
# -----------------------------------------------------------

## -- Start
message(" -- Start dashboard ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_start(${MODEL} TRACK ${MODEL})

## -- Update
message(" -- Update ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)

## -- Configure
message(" -- Configure ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_configure(BUILD  "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)

## -- Build
message(" -- Build ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_build(BUILD  "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)


# -----------------------------------------------------------
# -- Actual testing
# -----------------------------------------------------------
# The tests are in another repository and will run using those
# scripts. CTest reports on the test vouchers

## -- Clone the test repo if necessary
if(NOT EXISTS "${NALURTEST_DIR}")
       execute_process(COMMAND "${CTEST_GIT_COMMAND} clone ${NALURTEST_REPO_URL} ${NALURTEST_DIR}"
                       WORKING_DIRECTORY ${NIGHTLY_DIR} )
endif(NOT EXISTS "${NALURTEST_DIR}")

## -- Clean the test vouchers
file(GLOB old_vouchers "${NIGHTLY_TEST_RESULT_DIRECTORY}/*/PASS")
file(REMOVE ${old_vouchers})

## -- Run the tests using the script
message(" -- NaluRtest ${MODEL} - ${CTEST_BUILD_NAME} --")
execute_process(COMMAND "./run_tests.sh"
                WORKING_DIRECTORY ${NALURTEST_DIR}
                RESULT_VARIABLE nalurtestres
                OUTPUT_VARIABLE nalurtestlog
                ERROR_VARIABLE nalurtestlog)
file(WRITE ${NIGHTLY_DIR}/NaluRtest.txt "${nalurtestlog}")

## -- Run CTest 
message(" -- Test ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_test(BUILD  "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)

## -- SUBMIT
message(" -- Submit ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_submit( RETURN_VALUE res)

message(" -- Finished ${MODEL}  - ${CTEST_BUILD_NAME} --")
