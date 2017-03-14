Testing Nalu
============

Nalu's regression tests and unit tests are run nightly using the GCC and Intel 
compilers against the Trilinos master and development branches on a machine 
at NREL. The results can be seen at the `CDash Nalu website <http://my.cdash.org/index.php?project=Nalu>`__.


Run Tests Locally
-----------------

The nightly tests are implemented using `CTest <https://cmake.org/cmake/help/v3.7/manual/ctest.1.html>`__ and
are available to developers to run locally as well. Due to the nature of error propagation of 
calculations in computers, results of regression testing can be difficult to keep consistent. 
Output from Nalu can vary from established reference data for the regression tests based on the compiler you 
are using, the types of optimizations set, and the versions of the third-party libraries Nalu 
utilizes, along with the blas/lapack implementation in use. Therefore it may make sense when 
you first checkout Nalu, to create your own reference data for the tests for the machine and 
configuration you are using. The following instructions will describe how to run Nalu's tests.

Since Nalu's tests require a large amount of data, this data is hosted in a separate repository 
from Nalu. The first step is to checkout this repo. We will be checking out in our ``${HOME}`` directory.

::

   cd ${HOME} && git clone https://github.com/NaluCFD/NaluRtest.git

Once this repo is cloned, you will need to configure Nalu with testing on and with the location 
of this repo and another directory where the tests will be run. To configure Nalu with testing 
enabled, in Nalu's existing ``build`` directory, we will run this command:

::

   cmake -DTrilinos_DIR:PATH=`spack location -i nalu-trilinos` \
         -DYAML_DIR:PATH=`spack location -i yaml-cpp` \
         -DENABLE_INSTALL:BOOL=ON \
         -DCMAKE_BUILD_TYPE=RELEASE \
         -DENABLE_TESTS:BOOL=ON \
         -DNALURTEST_DIR:PATH=${HOME}/NaluRtest \
         -DRUNNALURTEST_DIR:PATH=${HOME}/runNaluRtest \
         ..

Note we have chosen to originally build Nalu with Spack in our case, hence the use 
of ``spack location -i <package>`` to locate our Yaml and Trilinos installations. 
So we use ``-DENABLE_TESTS:BOOL=ON`` to enable CTest, ``-DNALURTEST_DIR:PATH=${HOME}/NaluRtest`` 
to tell Nalu where the ``NaluRtest`` repo is, and finally ``-DRUNNALURTEST_DIR:PATH=${HOME}/runNaluRtest`` 
to tell Nalu where it is allowed to run the tests. The ``runNaluRtest`` directory will be created if it 
doesn't exist during testing, and this directory will also be erased each time the tests are run. Once 
Nalu is configured, you should be able to run the tests by building Nalu in the ``build`` directory, 
and running ``make test``.

There are advantages to using CTest here, such as being able to run subsets of the tests, or tests 
matching a particular regular expression for example. To do so, in the ``build`` directory, you can run 
``ctest -R femHC`` to run the test matching the ``femHC`` regular expression. See ``ctest -h`` for 
more options in running tests.

To define your own tolerance for tests, at configure time, add ``-DTEST_TOLERANCE=0.0001`` for example 
to the Nalu CMake configure line.

Updating Reference Data for Your Machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When running the tests, the norms for each test are the output and they are 'diffed' against 
the 'gold' norms that we have established for each test. To dictate whether or not test passes, 
we use a chosen tolerance in which we allow the results to deviate from the 'gold' norm.  As stated 
earlier, these 'gold' norms do not reflect every configuration of Nalu, per compiler, optimization, 
TPL versions, blas/lapack version, etc. This tolerance is currently defined in the ``CMakeLists.txt`` 
in Nalu's ``reg_tests`` directory. This tolerance can also be passed into Nalu at configure time using 
``-DTEST_TOLERANCE=0.0000001`` for example. To update the 'gold' norms to your configuration, merely 
run the tests once, and copy the ``*.norm`` files in ``runNaluRtest`` to the corresponding test location 
in ``NaluRtest`` to overwrite the current 'gold' norms.

Adding Tests to Nalu
--------------------

.. _add-test:

The testing infrastructure is almost completely confined to the ``reg_tests`` directory. To add a test 
to Nalu, we need to most likely modify two files. First we modify ``CTestPrepareTests.cmake.in`` to 
copy the files and setup the directory in ``runNaluRtest`` necessary for our new test. Make sure to copy 
the ``*.xml`` file if necessary. This CMake file is automatically run before the tests are run to prepare 
the directories they run in. It should be sufficient to use the same structure as the existing 
tests for the new test to accomplish this. Once the test has preparation steps established, the test 
itself can be added to the list of CTest tests by adding the test to the ``CTestList.cmake`` file. 
For a single regression test, provided it is uniform like the existing tests, it can likely be 
added by adding a single line with the test name and amount of processes you would like to run the 
test with. For example:

::

    add_test_r(mytest 6)

After this has been done, in the ``NaluRtest`` repo, you should add a directory corresponding to your 
test name and include the input file, ``mytest.i``, and reference output file ``mytest.norm.gold``.

To see commands used when running the tests, see the functions in the ``CTestList.cmake`` file. These 
functions ultimately create ``CTestTestFile.cmake`` files in the CMake build directory at configure time. 
You can see the exact commands used for each test after configure in the 
``build/reg_tests/CTestTestFile.cmake`` file.
