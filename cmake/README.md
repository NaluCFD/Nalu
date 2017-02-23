
# This is out-of-date at the moment

Testing is enabled through CTest and CDash and the results are posted
[here](http://my.cdash.org/index.php?project=Nalu).

## Prerequisites

- All third party libraries have been built (including Trilinos).

- The do-configNalu scripts point to the right directories. This
  should be changed at some point. Maybe SPACK can solve this.
    
## Running
	  
The CTest script can be run using the following command, for example: 
```{bash}
cd Nalu/build
ctest -DNIGHTLY_DIR=$nightlyDirectory -S ../cmake/ctest_peregrine_gcc.cmake
```
where `$nightlyDirectory` is the directory where the tests are run.

## Adding a new platform

If you want to add a new testing platform, create a new file called,
for example, `ctest_platform_compiler.cmake` in the `cmake`
directory. In that file, define the compiler that should be used and
the Nalu configure script for that platform as such:
```
set(compiler                            "gcc")
set(CONFIG_FILE                         "do-configNalu_script")
include(../cmake/ctest_common.cmake)
```
You can look at `ctest_peregrine_gcc.cmake` for inspiration.


## Adding a new test to the CTest suite

There are two files that must be edited when adding a test to the
regression test suite, both located in the `cmake` directory:

- `ctest_prepare_tests.cmake`: sets up the directory structure needed for the test;

- `CTestTestfile.cmake`: defines the test and its properties.

This presumes you have created a directory for your test in the
NaluRtest repo with the appropriate input files, etc.

### `ctest_prepare_tests.cmake`

This file sets up the directory structure needed for the test. It
copies files from the NaluRtest directory into a local testing
directory. You are responsible for ensuring that the proper files are
copied over.

A simple entry in this file would be something like:
```
file(COPY               "${NALURTEST_NIGHTLY_DIR}/periodic3dElem"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_MESH_DIR}/periodic3d.g"
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/periodic3dElem")
```

This set of commands copies the test directory to a local nightly
directory (including all the files inside of it) and it copies some
extra files located in other directories (such as the mesh and
`milestone.xml`).


A more complicated set of instructions (with file globbing) would look like:
```
file(COPY               "${NALURTEST_NIGHTLY_DIR}/dgNonConformal3dFluidsHexTet"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/twoBlockMeshHexTet_cgs.g.*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformal3dFluidsHexTet")
```

### `CTestTestfile.cmake`

This file defines the test and its properties.

A sample entry would look like
```
add_test(periodic3dElemNp4 "periodic3dElemNp4.sh")
set_tests_properties(periodic3dElemNp4
   PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
               FAIL_REGULAR_EXPRESSION "FAILED"
   TIMEOUT 400
   PROCESSORS 4
   WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/periodic3dElem")
```

This defines, in order:

- the test exectutable file (e.g. `periodic3dElemNp4.sh`)

- the regular expression to look for the test status

- the maximum time in seconds run the test before sending an error

- the number of processors to use

- the directory where the test should be run
