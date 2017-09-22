Building Nalu Manually
======================

If you prefer not to build using Spack, below are instructions which
describe the process of building Nalu by hand.

Linux
-----

The Linux platform installation instructions follow. Compilers and MPI
are expected to be already installed. If they are not, please follow the
open-mpi build instructions. Below, we are using openmpi-1.8.8 and
gcc-4.8.5

Setup
~~~~~

Prepare the third-party library (TPL) build process by defining some 
code locations, e.g., ``gitHubWork/scratch_build`` and set some paths. 
One might choose to keep a ``nalu_module_4.7.2`` file.

::

    mkdir <your_base_dir>
    cd <your_base_dir>
    export nalu_build_dir=$PWD
    mkdir $nalu_build_dir/packages
    mkdir $nalu_build_dir/install
    mkdir $nalu_build_dir/install/lib
    export LD_LIBRARY_PATH=/YOUR_PATH_TO_MPI_LIB/1.6.4-gcc-4.7.2-RHEL6/lib:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=/YOUR_PATH_TO_GCC/gcc/4.7.2-RHEL6/lib64:$LD_LIBRARY_PATH
    PATH=/YOUR_PATH_TO_MPI/bin:$PATH
    PATH=$nalu_build_dir/install/bin:$PATH
    PATH=$nalu_build_dir/install/trilinos/bin:$PATH

CMake, Version 3.1.0
~~~~~~~~~~~~~~~~~~~~

CMake is provided by the `CMake team <http://www.cmake.org/download/>`__.

::

    cd $nalu_build_dir/packages
    curl -o cmake-3.1.0-rc2.tar.gz http://www.cmake.org/files/v3.1/cmake-3.1.0-rc2.tar.gz
    tar -xvf cmake-3.1.0-rc2.tar.gz

Build (we assume the configure command will find the needed compilers):

::

    cd $nalu_build_dir/packages/cmake-3.1.0-rc2
    ./configure --prefix=$nalu_build_dir/install
    gmake
    gmake install

SuperLU, Version 4.3
~~~~~~~~~~~~~~~~~~~~

SuperLU is provided by `here <http://crd-legacy.lbl.gov/~xiaoye/SuperLU/>`__.

::

    cd $nalu_build_dir/packages
    curl -o superlu_4.3.tar.gz http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_4.3.tar.gz
    tar -xvf superlu_4.3.tar.gz

Build:

::

    cd $nalu_build_dir/packages/SuperLU_4.3
    cp MAKE_INC/make.linux make.inc

To find out what the correct platform extension PLAT is:

::

    uname -m

Edit ``make.inc`` as shown below (diffs shown from baseline).

::

    PLAT = _x86_64
    SuperLUroot   = /your_path/install/SuperLU_4.3 i.e., $nalu_build_dir/install/SuperLU_4.3
    BLASLIB       = -L/usr/lib64 -lblas
    CC            = mpicc
    FORTRAN       = mpif77

On some platforms, the ``$nalu_build_dir`` may be mangled. In such
cases, you may need to use the entire path to ``install/SuperLU_4.3``.

Next, make some new directories:

::

    mkdir $nalu_build_dir/install/SuperLU_4.3
    mkdir $nalu_build_dir/install/SuperLU_4.3/lib
    mkdir $nalu_build_dir/install/SuperLU_4.3/include
    cd $nalu_build_dir/packages/SuperLU_4.3
    make
    cp SRC/*.h $nalu_build_dir/install/SuperLU_4.3/include

Libxml2, Version 2.9.2
~~~~~~~~~~~~~~~~~~~~~~

Libxml2 is found `here <http://www.xmlsoft.org/sources/>`__.

::

    cd $nalu_build_dir/packages
    curl -o libxml2-2.9.2.tar.gz http://www.xmlsoft.org/sources/libxml2-2.9.2.tar.gz
    tar -xvf libxml2-2.9.2.tar.gz

Build (note that python is not required and, hence, the ``-without-python``
option to config):

::

    cd $nalu_build_dir/packages/libxml2-2.9.2
    CC=mpicc CXX=mpicxx ./configure -without-python --prefix=$nalu_build_dir/install
    make
    make -k install

Boost, Version 1.55.0
~~~~~~~~~~~~~~~~~~~~~

Boost is found `here <http://www.boost.org>`__.

Directions for 1.55.0 are as follows:

::

    cd $nalu_build_dir/packages
    curl -o boost_1_55_0.tar.gz http://iweb.dl.sourceforge.net/project/boost/boost/1.55.0/boost_1_55_0.tar.gz
    tar -zxvf boost_1_55_0.tar.gz

Build:

::

    cd $nalu_build_dir/packages/boost_1_55_0

Note: There must be a space before the semicolon at the end of the
"using mpi" line in the ``user-config.jam`` file.

::

    echo "using mpi : `which mpicxx` ;" >> ./tools/build/v2/user-config.jam 
    ./bootstrap.sh --prefix=$nalu_build_dir/install --with-libraries=signals,regex,filesystem,system,mpi,serialization,thread,program_options,exception 
    ./b2 -j 4 2>&1 | tee boost_build_one
    ./b2 -j 4 install 2>&1 | tee boost_build_intall

For newer versions, e.g., 1.60.0:

::

        ./bootstrap.sh --prefix=$nalu_build_dir/install --with-libraries=signals,regex,filesystem,system,mpi,serialization,thread,program_options,exception

Next, edit ``project-config.jam`` and add a 'using mpi', e.g,

using mpi: /path/to/mpi/openmpi/bin/mpicc

::

        ./b2 -j 4 2>&1 | tee boost_build_one
        ./b2 -j 4 install 2>&1 | tee boost_build_intall

Yaml-cpp
~~~~~~~~

For versions of Nalu after the v1.1.0-release, Yaml is provided under
`Github <https://github.com/jbeder/yaml-cpp>`__.

::

    cd $nalu_build_dir/packages
    git clone https://github.com/jbeder/yaml-cpp

Build Yaml-cpp:

::

    cd $nalu_build_dir/packages/yaml-cpp
    mkdir build
    cd build
    cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CXX_FLAGS=-std=c++11 -DCMAKE_CC_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX=$nalu_build_dir/install ..
    make
    make install

Pre-v1.1.0-release; Yaml-cpp, Version 0.3.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For all versions of Nalu at, or before the v1.1.0-release, the formal
version of YAML is 0.3.0. There is no backward compatibility between the
versions of YAML.

Yaml is provided by
`code.google.com <https://code.google.com/p/yaml-cpp/downloads/detail?name=yaml-cpp-0.3.0.tar.gz&can=2&q=>`__.

::

    cd $nalu_build_dir/packages
    curl -o yaml-cpp-0.3.0.tar.gz https://yaml-cpp.googlecode.com/files/yaml-cpp-0.3.0.tar.gz
    tar -zxvf yaml-cpp-0.3.0.tar.gz
    mv yaml-cpp yaml-cpp-0.3.0

Build yaml-cpp:

::

    cd $nalu_build_dir/packages/yaml-cpp-0.3.0
    mkdir build
    cd build
    cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CC_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX=$nalu_build_dir/install ..
    make
    make install

Zlib, Version 1.2.8
~~~~~~~~~~~~~~~~~~~

Zlib is provided `here <http://www.zlib.net>`__.

::

    cd $nalu_build_dir/packages
    curl -o zlib-1.2.8.tar.gz http://zlib.net/zlib-1.2.8.tar.gz
    tar -zxvf zlib-1.2.8.tar.gz

Build Zlib:

::

    cd $nalu_build_dir/packages/zlib-1.2.8
    CC=gcc CXX=g++ CFLAGS=-O3 CXXFLAGS=-O3 ./configure --prefix=$nalu_build_dir/install/
    make
    make install

Hdf5, Version 1.8.12
~~~~~~~~~~~~~~~~~~~~

Hdf5 1.8.12 is provided by the
`HDF <http://www.hdfgroup.org/downloads/index.html>`__ group.

::

    cd $nalu_build_dir/packages/
    curl -o hdf5-1.8.12.tar.gz http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.12/src/hdf5-1.8.12.tar.gz
    tar -zxvf hdf5-1.8.12.tar.gz

Build (parallel enabled):

::

    cd $nalu_build_dir/packages/hdf5-1.8.12
    ./configure CC=mpicc FC=mpif90 CXX=mpicxx CXXFLAGS="-fPIC -O3" CFLAGS="-fPIC -O3" FCFLAGS="-fPIC -O3" --enable-parallel --with-zlib=$nalu_build_dir/install --prefix=$nalu_build_dir/install
    make
    make install
    make check
        

Full Parallel-Enabled Nalu using NetCDF (V. 4.3.3.1) and Parallel NetCDF (V. 1.6.1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to support all aspects of Nalu's parallel models, this
combination of products is required.

Parallel NetCDF, Version 1.6.1
******************************

Parallel NetCDF is provided on the `Argon Trac
Page <https://trac.mcs.anl.gov/projects/parallel-netcdf/wiki/Download>`__.

::

    cd $nalu_build_dir/packages/
    tar -zxvf parallel-netcdf-1.6.1.tar.gz

Configure, build and install:

::

    cd parallel-netcdf-1.6.1
    ./configure --prefix=$nalu_install_dir CC=mpicc FC=mpif90 CXX=mpicxx CFLAGS="-I$nalu_install_dir/include -O3" LDFLAGS=-L$nalu_install_dir/lib --disable-fortran
    make
    make install

Note that we have created an install directory that might look like ``$nalu_build_dir/install``.

NetCDF Version 4.3.3.1
**********************

NetCDF is provided on
`Github <https://github.com/Unidata/netcdf-c/releases>`__.

::

    cd $nalu_build_dir/packages/
    curl -o netcdf-c-4.3.3.1.tar.gz https://codeload.github.com/Unidata/netcdf-c/tar.gz/v4.3.3.1
    tar -zxvf netcdf-c-4.3.3.1.tar.gz

Configure, build and install:

::

    cd netcdf-c-4.3.3.1
    ./configure --prefix=$nalu_install_dir CC=mpicc FC=mpif90 CXX=mpicxx CFLAGS="-I$nalu_install_dir/include -O3" LDFLAGS=-L$nalu_install_dir/lib --enable-pnetcdf --enable-parallel-tests --enable-netcdf-4 --disable-shared --disable-fsync --disable-cdmremote --disable-dap --disable-doxygen --disable-v2
    make -j 4 
    make install
    make check

Note that when using Parallel NetCDF, the proper install directories
must be added to the Trilinos configuration file.

Partial Parallel-Enabled Nalu using NetCDF, Version 4.3.1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If Parallel NetCDF is omitted, follow the instructions below. On some
platforms, however, autodecompostion may fail.

NetCDF is provided on
`Github <https://github.com/Unidata/netcdf-c/releases>`__.

Scroll down until you see "NetCDF-C 4.3.1.1 (Bugfix Release)" or similar
Click on the "Source (tar.gz)" button to download and then move the tar
file to:

::

    cd $nalu_build_dir/packages/
    curl -o netcdf-c-4.3.1.1.tar.gz https://codeload.github.com/Unidata/netcdf-c/tar.gz/v4.3.1.1
    tar -zxvf netcdf-c-4.3.1.1.tar.gz

Possibly, 4.3.1.1 is hard to get... If so, use the following:

::

    curl -o netcdf-c-4.3.1-rc2.tar.gz https://codeload.github.com/Unidata/netcdf-c/tar.gz/v4.3.1-rc2

Complex Models (expert usage only)
**********************************

In ``netcdf/include/netcdf.h``, the following defines need to be changed to
support complex models.

::

    #define NC_MAX_DIMS     65536    /* max dimensions per file */
    #define NC_MAX_VARS     524288   /* max variables per file */

For a definiton of Complex Models, please note the following page:

`complexModels <https://github.com/gsjaardema/seacas/blob/master/NetCDF-Mapping.md>`__

Care should be taken with these settings as sometimes the above setting
can exceed platform resources and, therefore, casue fails in the
installation test suite.

Build (with parallel I/O):

::

    cd $nalu_build_dir/packages/netcdf-c-4.3.1.1
    ./configure --prefix=$nalu_build_dir/install CC=mpicc FC=mpif90 CXX=mpicxx CFLAGS="-I$nalu_build_dir/install/include -O3" LDFLAGS=-L$nalu_build_dir/install/lib --disable-fsync --disable-cdmremote --disable-dap --disable-shared --disable-doxygen
    make -j 4 
    make install
    make check

Trilinos
~~~~~~~~

Trilinos is managed by the `Trilinos <http://www.trilinos.org>`__
project and can be found on Github.

Clone the latest version of Trilinos within
``$nalu_build_dir/packages``:

::

    cd $nalu_build_dir/packages/
    git clone https://github.com/trilinos/Trilinos.git

In some cases, the master Trilinos code base may have build issues. This
is a rare occurance, however, some aspects to Trilinos that Nalu
require, e.g., Tpetra, kokkos, STK and Muelu are in *active*
development. If problems arise, one can revert back to a possible
successful SHA-1 using bisect. Again, this is hopefully going to be
mitigated by the strong SQA efforts at SNL.

Nalu Releases
*************

Unfortunately, Github does not allow for a "live" wiki for each of the
existing branches of Nalu.wiki. As such, instructions for the particular
releases have been embedded within this head wiki file.

Release v1.0.0-release
^^^^^^^^^^^^^^^^^^^^^^

For the formal Nalu v1.0.0-release, checkout the following Trilinos
Version:

::

        git checkout trilinos-release-12-0-branch   

This version is the expected Trilinos code base for the v1.0.0-release
Nalu code base. Now proceed to the build section.

Head Code Base
^^^^^^^^^^^^^^

Proceed to the build section without checking out the Trilinos
12-0-branch.

Build
*****

Create new folder in Trilinos called build:

::

    cd $nalu_build_dir/packages/Trilinos
    mkdir build

Place into build the script one of the ``do-configTrilinos_*`` files.

``do-configTrilinos_*`` will be used to run cmake to build trilinos
correctly for Nalu. Note that there are two files: one for 'release'
and the other 'debug'. The files can be found on the Nalu GitHub site
or copied from ``$nalu_build_dir/packages/Nalu/build``, which is
created in the Nalu build step documented below. For example:

Pull latest version of ``do-configTrilinos_*`` from Nalu's GitHub site:

::

    curl -o $nalu_build_dir/packages/Trilinos/build/do-configTrilinos_release https://raw.githubusercontent.com/NaluCFD/Nalu/master/build/do-configTrilinos_release

Or if you create the Nalu directory as directed below, simply copy one
of the ``do-configTrilinos_*`` files from local copy of Nalu's git
repository:

::

    cp $nalu_build_dir/packages/Nalu/build/do-configTrilinos_release $nalu_build_dir/packages/Trilinos/build

Now edit ``do-configTrilinos_release`` to modify the paths so they point
to ``$nalu_build_dir/install``.

::

    cd $nalu_build_dir/packages/Trilinos/build
    chmod +x do-configTrilinos_release

Make sure all other paths to netcdf, hdf5, etc., are correct (in
addition to open-mpi).

::

    ./do-configTrilinos_release
    make
    make install

If after the make, one notes issues with hdf5 and netcdf references not
found, add the following:

::

    -DTPL_Netcdf_LIBRARIES:PATH="${netcdf_install_dir}/lib/libnetcdf.a;${hdf_install_dir}/lib/libhdf5_hl.a;${hdf_install_dir}/lib/libhdf5.a;${z_install_dir}/lib/libz.a" \

Just below the netcdf option within the Seacas do-config sections:

::

    -DTPL_ENABLE_Netcdf:STRING=ON \ 

ParaView Catalyst
~~~~~~~~~~~~~~~~~

Optionally enable `ParaView Catalyst <https://www.paraview.org/in-situ/>`__
for in-situ visualization with Nalu. These instructions can be skipped if 
you do not require in-situ visualization with Nalu.

Build ParaView SuperBuild v5.3.0
********************************

The `ParaView SuperBuild <https://gitlab.kitware.com/paraview/paraview-superbuild>`__ 
builds ParaView along with all dependencies necessary to enable Catalyst with Nalu.
Clone the ParaView SuperBuild within ``$nalu_build_dir/packages``:

::

    cd $nalu_build_dir/packages/
    git clone --recursive https://gitlab.kitware.com/paraview/paraview-superbuild.git
    cd paraview-superbuild
    git fetch origin
    git checkout v5.3.0
    git submodule update

Create a new build folder in ``$nalu_build_dir/``:

::

    cd $nalu_build_dir
    mkdir paraview-superbuild-build
    cd paraview-superbuild-build

Copy ``do-configParaViewSuperBuild`` to ``paraview-superbuild-build``.
Edit ``do-configParaViewSuperBuild`` to modify the defined paths as
follows:

::

    mpi_base_dir=<same MPI base directory used to build Trilinos>
    nalu_build_dir=<path to root nalu build dir>

Make sure the MPI library names are correct.

::

    ./do-configParaViewSuperBuild
    make -j 8
   
Build Nalu ParaView Catalyst Adapter
************************************

Create a new build folder in ``$nalu_build_dir/``:

::

    cd $nalu_build_dir
    mkdir nalu-catalyst-adapter-build
    cd nalu-catalyst-adapter-build

Copy ``do-configNaluCatalystAdapter`` to ``nalu-catalyst-adapter-build``.
Edit ``do-configNaluCatalystAdapter`` and modify ``nalu_build_dir`` at the
top of the file to the root build directory path.

::

    ./do-configNaluCatalystAdapter
    make
    make install

Nalu
~~~~

Nalu is provided `here <https://github.com/NaluCFD/Nalu>`__.

No doubt, you already have cloned Nalu. If not, execute the following
command in the location that you want Nalu:

::

    git clone https://github.com/NaluCFD/Nalu.git

Nalu Releases
*************

One may either build the released Nalu version, v1.0.0-release, or the
head code base.

Release v1.0.0-release
^^^^^^^^^^^^^^^^^^^^^^

For the formal Nalu v1.0.0-release, you should have already cloned
Trilinos and built the 12.0 release version of Trilinos. To obtain the
consistent Nalu version, after the clone, checkout the Nalu release:

::

    git checkout v1.0.0-release

Now proceed to the build section below.

Head Code Base
^^^^^^^^^^^^^^

Proceed to the build section without checking out the Nalu
v1.0.0-release code repository.

Build
*****

In ``Nalu/build``, you will find the
`CMakeLists.txt <https://github.com/NaluCFD/Nalu/blob/master/CMakeLists.txt>`__
and
`do-configNalu <https://github.com/NaluCFD/Nalu/blob/master/build/do-configNalu_release>`__.

Copy the ``do-configNalu_release`` or ``debug`` file to a new, non-tracked
file:

::

    cp do-configNalu_release do-configNaluNonTracked

Edit the paths at the top of the files by defining the
``nalu_build_dir variable``. Within ``Nalu/build``, execute the
following commands:

::

    ./do-configNaluNonTracked
    make 

This process will create ``naluX`` within the ``Nalu/build`` location.
One may need to create a new ``yaml-cpp`` directory (and copy ``src/include``
files). You may also build a debug executable by modifying the Nalu
config file to use "Debug". In this case, a ``naluXd`` executable is
created.

Other useful tools from, e.g., seacas, are under
``/usr/local/packages/install/trilinos/bin``

Build Nalu with ParaView Catalyst Enabled
*****************************************

If you have built ParaView Catalyst and the Nalu ParaView Catalyst Adapter, you
can build Nalu with Catalyst enabled.

In ``Nalu/build``, find ``do-configNaluCatalyst``. Copy ``do-configNaluCatalyst`` to
a new, non-tracked file:

::

    cp do-configNaluCatalyst do-configNaluCatalystNonTracked
    ./do-configNaluCatalystNonTracked
    make 

The build will create the same executables as a regular Nalu build, and will also create a  
bash shell script named ``naluXCatalyst``.  Use ``naluXCatalyst`` to run Nalu
with Catalyst enabled.  It is also possible to run ``naluX`` with Catalyst enabled by
first setting the environment variable:

::

   export CATALYST_ADAPTER_INSTALL_DIR=$nalu_build_dir/install

Nalu will render images to Catalyst in-situ if it encounters the keyword ``catalyst_file_name``
in the ``output`` section of the Nalu input deck. The ``catalyst_file_name`` command specifies the
path to a text file containing ParaView Catalyst input deck commands. Consult the ``catalyst.txt`` files
in the following Nalu regression test directories for examples of the Catalyst input deck command syntax:

::

    ablForcingEdge/
    mixedTetPipe/
    steadyTaylorVortex/

::

    output:
      output_data_base_name: mixedTetPipe.e
      catalyst_file_name: catalyst.txt

When the above regression tests are run, Catalyst is run as part of the regression test. The regression
test checks that the correct number of image output files have been created by the test.

The Nalu Catalyst integration also supports running Catalyst Python script files exported from the ParaView GUI.
The procedure for exporting Catalyst Python scripts from ParaView is documented in the 
`Catalyst user guide <https://www.paraview.org/in-situ/>`__. To use an exported Catalyst script, insert 
the ``paraview_script_name`` keyword in the ``output`` section of the Nalu input deck. The argument for
the ``paraview_script_name`` command contains a file path to the exported script. 

::

    output:
      output_data_base_name: mixedTetPipe.e
      paraview_script_name: paraview_exported_catalyst_script.py


Mac OS X
--------

The Mac platform installation is managed by the extensive usage of
"Homebrew". This package provides many Mac OS X builds and installations
required for Nalu.

Homebrew
~~~~~~~~

Download and Install
`Homebrew <https://github.com/Homebrew/homebrew/wiki/Installation>`__ on
your local home terminal:

Now you can install Homebrew; ``brew doctor`` will be the first line
command required.

::

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    brew doctor

Packages to be obtained from Homebrew
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OpenMpi
*******

::

    brew install openmpi

Cmake
*****

::

    brew search cmake  //cmake is there and will pop up on the next line
    brew install cmake

Libxml2
*******

::

    brew search libxml2
    brew install libxml2

Boost, 1.60.0_2
***************

::

    brew install boost

SuperLU, 4.3
************

::

    brew tap homebrew/science
    brew search superlu 
    brew install superlu43

The latest version of SuperLU provided by Homebrew is not compatible
with head Trilinos due to Trilinos' usage of some deprecated methods.

Once done with the Homebrew install, make sure everybody can read the files:

::

    sudo chmod -R u+rwX,go+rX,go-w /usr/local/Cellar

Non-Homebrew
~~~~~~~~~~~~

Other Nalu required libraries must be managed outside of the Homebrew
environment.

Helpful Notes
*************

How to specify install location and compilers with ``configure`` and
``cmake``:

When building packages that use configure, do:

::

    ./configure PREFIX=/usr/local/packages/install CC=mpicc CXX=mpicxx

When building packages that use cmake, do a ``mkdir build; cd build``
then

::

    cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CC_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX:PATH=/myPath/install ..`

For all non-Homebrew packages, ``myPath`` above will be ``/usr/local/packages``.

::

    mkdir /usr/local/packages
    mkdir /usr/local/packages/install

Yaml-cpp
********

For versions of Nalu after the v1.1.0-release, Yaml is provided under
`Github <https://github.com/jbeder/yaml-cpp>`__

::

    cd $nalu_build_dir/packages
    git clone https://github.com/jbeder/yaml-cpp

Build Yaml-cpp:

::

    cd $nalu_build_dir/packages/yaml-cpp
    mkdir build
    cd build
    cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CXX_FLAGS=-std=c++11 -DCMAKE_CC_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX=$nalu_build_dir/install ..
    make
    make install

Pre-v1.1.0-release; Yaml-cpp, Version 0.3.0
*******************************************

Yaml is provided by
`code.google.com <https://code.google.com/p/yaml-cpp/downloads/detail?name=yaml-cpp-0.3.0.tar.gz&can=2&q=>`__.

Follow the yaml installation directions. This process will put the
created files in ``/user/local/include/yaml-cpp``. Below are some high level
points:

::

    mkdir /usr/local/packages
    cd packages/
    curl -o yaml-cpp-0.3.0.tar.gz https://yaml-cpp.googlecode.com/files/yaml-cpp-0.3.0.tar.gz 
    tar -zxvf yaml-cpp-0.3.0.tar.gz 
    mv yaml-cpp yaml-cpp-0.3.0

This series of commands will create ``/usr/local/packages/yaml-cpp-0.3.0``.

Next, build Yaml-cpp:

::

    cd /usr/local/packages/yaml-cpp-0.3.0
    mkdir build
    cd build
    cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CC_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX:PATH=/usr/local/packages/install ..
    make
    make install

Zlib, 1.2.8
***********

Zlib is provided by the `zlib <http://www.zlib.net>`__ project.

::

    cd /usr/local/packages/
    curl -o zlib-1.2.8.tar.gz http://zlib.net/zlib-1.2.8.tar.gz
    tar -zxvf zlib-1.2.8.tar.gz 

Build Zlib:

::

    cd /usr/local/packages/zlib-1.2.8
    CC=gcc CXX=g++ CFLAGS=-O3 CXXFLAGS=-O3 ./configure --archs="-arch x86_64" --prefix=/usr/local/packages/install/
    make
    make install

Hdf5, 1.8.12
************

Hdf5 1.8.12 is provided by the
`HDF <http://www.hdfgroup.org/downloads/index.html>`__ group.

::

    cd /usr/local/packages/
    curl -o hdf5-1.8.12.tar.gz http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.12/src/hdf5-1.8.12.tar.gz
    tar -zxvf hdf5-1.8.12.tar.gz 

This series of commands will create ``/usr/local/packages/hdf5-1.8.12``.

Build:

::

    cd /usr/local/packages/hdf5-1.8.12
    ./configure CC=mpicc FC=mpif90 CXX=mpicxx CXXFLAGS="-fPIC -O3" CFLAGS="-fPIC -O3" FCFLAGS="-fPIC -O3" --enable-parallel --with-zlib=/usr/local/packages/install --prefix=/usr/local/packages/install
    make
    make install
    make check

Full Parallel-Enabled Nalu using NetCDF (V. 4.3.3.1) and Parallel NetCDF (V. 1.6.1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to support all aspects of Nalu's parallel models, this
combination of products is required.

Parallel NetCDF, Version 1.6.1
******************************

Parallel NetCDF is provided on the `Argon Trac
Page <https://trac.mcs.anl.gov/projects/parallel-netcdf/wiki/Download>`__.

::

    cd $nalu_build_dir/packages/
    tar -zxvf parallel-netcdf-1.6.1.tar.gz

Configure, build and install:

::

    cd parallel-netcdf-1.6.1
    ./configure --prefix=/usr/local/packages/install CC=mpicc FC=mpif90 CXX=mpicxx CFLAGS="-I/usr/local/packages/install/include -O3" LDFLAGS=-L/usr/local/packages/install/lib --disable-fortran
    make
    make install

Note that we have created an install directory that might look like:
``$nalu_build_dir/install``.

NetCDF Version 4.3.3.1
**********************

NetCDF is provided on
`Github <https://github.com/Unidata/netcdf-c/releases>`__.

::

    cd $nalu_build_dir/packages/
    curl -o netcdf-c-4.3.3.1.tar.gz https://codeload.github.com/Unidata/netcdf-c/tar.gz/v4.3.3.1
    tar -zxvf netcdf-c-4.3.3.1.tar.gz

Configure, build and install:

::

    cd netcdf-c-4.3.3.1
    ./configure --prefix=$nalu_install_dir CC=mpicc FC=mpif90 CXX=mpicxx CFLAGS="-I/usr/local/packages/install/include -O3" LDFLAGS=-L/usr/local/packages/install/lib --enable-pnetcdf --enable-parallel-tests --enable-netcdf-4 --disable-shared --disable-fsync --disable-cdmremote --disable-dap --disable-doxygen --disable-v2
    make -j 4 
    make install
    make check

Note that when using Parallel NetCDF, the proper install directories
must be added to the Trilinos configuration file.

Partial Parallel-Enabled Nalu using NetCDF, Version 4.3.1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If Parallel NetCDF is omitted, follow the instructions below. On some
platforms, however, autodecompostion may fail.

Netcdf is provided on
`Github <https://github.com/Unidata/netcdf-c/releases>`__.

Scroll down until you see "NetCDF-C 4.3.1.1 (Bugfix Release)" or similar
Click on the "Source (tar.gz)" button to download and then move the tar
file to: ``/usr/local/packages``

::

    cd /usr/local/packages
    tar -xvf netcdf-c-4.3.1.1.tar

This process will create ``/usr/local/packages/netcdf-c-4.3.1.1``.

Build:

::

    cd /usr/local/packages/netcdf-c-4.3.1.1
    ./configure --prefix=/usr/local/packages/install CC=mpicc FC=mpif90 CXX=mpicxx CFLAGS="-I/usr/local/packages/install/include -O3" LDFLAGS="-L/usr/local/packages/install/lib" --disable-fsync --disable-cdmremote --disable-dap --disable-shared --disable-doxygen
    make
    make install
    make check

Trilinos
********

Trilinos is managed by the `Trilinos <http://www.trilinos.org>`__
project and can be found on github.

Clone the latest version of Trilinos within ``/packages``:

::

    cd /usr/local/packages/
    git clone https://github.com/trilinos/Trilinos.git

In some cases, the master Trilinos code base may have build issues. This
is a rare occurance, however, some aspects to Trilinos that Nalu
require, e.g., Tpetra, kokkos, STK and Muelu are in *active*
development. If problems arise, one can revert back to a possible
successful SHA-1 using bisect. Again, this is hopefully going to be
mitigated by the strong SQA efforts at SNL.

Nalu Releases
~~~~~~~~~~~~~

Nalu has two main release so far.

Release v1.0.0-release
**********************

For the formal v1.0.0-release, check-out the following Trilinos Version:

::

        git checkout trilinos-release-12-0-branch       

This version is the expected Trilinos code base for the v1.0.0-release
Nalu code base. Now proceed to the build section.

Head Code Base
**************

Proceed to the build section without checking out the Trilinos
12-0-branch.

Build
*****

Create new folder in Trilinos called ``build``:

::

    cd /usr/local/packages/Trilinos
    mkdir build

Place into build the script one of the ``do-configTrilinos_*`` files.

``do-configTrilinos_*`` will be used to run cmake to build trilinos
correctly for Nalu. Note that there are two files: one for 'release'
and the other 'debug'. The files can be found on the Nalu GitHub site
or copied from ``$nalu_build_dir/packages/Nalu/build``, which is
created in the Nalu build step documented below. For example:

Pull latest version of ``do-configTrilinos_*`` from Nalu's GitHub site:

::

    curl -o $nalu_build_dir/packages/Trilinos/build/do-configTrilinos_release https://raw.githubusercontent.com/NaluCFD/Nalu/master/build/do-configTrilinos_release

Or if you create the Nalu directory as directed below, simply copy one
of the ``do-configTrilinos_*`` files from local copy of Nalu's git
repository:

::

    cp $nalu_build_dir/packages/Nalu/build/do-configTrilinos_release $nalu_build_dir/packages/Trilinos/build

Now edit the ``do-configTrilinos_release`` to modify the defined paths as
follows:

::

    mpi_base_dir=/usr/local/Cellar/open-mpi/1.8.3
    nalu_build_dir=/usr/local/packages

Next, note that some packages, i.e., boost and superLu were provided by
Homebrew. As such, make sure that ``boost_dir`` and ``super_lu`` point to
``/usr/local/Cellar``.

Build:

::

    cd /usr/local/packages/Trilinos/build
    ./do-configTrilinos_release
    make 
    make install

Nalu
~~~~

Nalu is provided `here <https://github.com/NaluCFD/Nalu>`__.

No doubt, you already have cloned Nalu. If not, execute the following
command in the location that you want Nalu:

::

    git clone https://github.com/NaluCFD/Nalu.git

Nalu Releases
*************

One may either build the released Nalu version, v1.0.0-release, or the
head code base.

Release v1.0.0-release
**********************

For the formal Nalu v1.0.0-release, you should have already cloned
Trilinos and built the 12.0 release version of Trilinos. To obtain the
consistent Nalu version, after the clone, checkout the Nalu release,

::

    git checkout v1.0.0-release

Now proceed to the build section below.

Head Code Base
**************

Proceed to the build section without checking out the Nalu
v1.0.0-release code repository.

Build
*****

In ``Nalu/build``, you will find the
`CMakeLists.txt <https://github.com/NaluCFD/Nalu/blob/master/CMakeLists.txt>`__
and
`do-configNalu_release <https://github.com/NaluCFD/Nalu/blob/master/build/do-configNalu_release>`__.

Again, note that there is a debug version as well. Copy the
``do-configNalu_release`` to a new, non-tracked file,

::

    cp do-configNalu_release do-configNaluNonTracked

Edit the paths at the top of the files by defining the
``nalu_build_dir variable`` as:

::

    nalu_build_dir=/usr/local/packages

Within ``Nalu/build``, execute the following commands:

::

    ./do-configNaluNonTracked
    make 

This process will create ``naluX`` within the ``Nalu/build`` location.
Setting the DEBUG CMake option will create a naluXd executable.

Other useful tools from, e.g., seacas, are under
``/usr/local/packages/install/trilinos/bin``.

