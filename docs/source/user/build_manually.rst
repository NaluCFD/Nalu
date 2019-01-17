Building Nalu Manually
======================

If you prefer not to build using Spack, below are instructions which describe the process of building Nalu by hand.

Linux and OSX
-------------

The instructions for Linux and OSX are mostly the same, except on each OS you may be able to use a package manager to 
install some dependencies for you. Using Homebrew on OSX is one option listed below. Compilers and MPI are expected to 
be already installed. If they are not, please follow the open-mpi build instructions. Start by creating a ``$nalu_build_dir`` 
and ``$nalu_install_dir`` in which to work. The former is the location in which the TPL source code resides while the 
later is for installation. A sample install might be ``/home/nalu_user/gitHubWork/scratch_build/install/gcc7.2.0`` while 
the build location, ``/home/nalu_user/gitHubWork/scratch_build``.


Homebrew
~~~~~~~~

If using OSX, you can install many dependencies using Homebrew. Install `Homebrew <https://github.com/Homebrew/homebrew/wiki/Installation>`__ 
on your local machine and reference the list below for some packages Homebrew can install for you which allows you to skip 
the steps describing the build process for each application, but not that you will need to find the location of the applications 
in which Homebrew has installed them, to use when building Trilinos and Nalu.

::

    brew install openmpi
    brew install cmake
    brew install libxml2
    brew install boost
    brew tap homebrew/science
    brew install superlu43


CMake v3.12.3
~~~~~~~~~~~~

CMake is provided `here <http://www.cmake.org/download/>`__.

Prepare:

::

    cd $nalu_build_dir/packages
    tar xf cmake-3.12.3.tar.gz

Build:

::

    cd $nalu_build_dir/packages/cmake-3.12.3
    ./configure --prefix=$nalu_install_dir/cmake/3.12.3
    make
    make install

SuperLU v4.3
~~~~~~~~~~~~

SuperLU is provided `here <http://crd-legacy.lbl.gov/~xiaoye/SuperLU/>`__.

Prepare:

::

    cd $nalu_build_dir/packages
    curl -o superlu_4.3.tar.gz http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_4.3.tar.gz
    tar xf superlu_4.3.tar.gz

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
    SuperLUroot   = /your_path_to_install/SuperLU_4.3 i.e., $nalu_install_dir/SuperLU/4.3
    BLASLIB       = -L/usr/lib64 -lblas
    CC            = mpicc
    FORTRAN       = mpif77

On some platforms, the ``$nalu_insall_dir`` may be mangled and, thus the make will fail. In such cases, you 
need to use the entire path to ``your_path_to_install/SuperLU_4.3``.

Next, make some new directories:

::

    mkdir $nalu_install_dir/SuperLU/4.3
    mkdir $nalu_install_dir/SuperLU/4.3/lib
    mkdir $nalu_install_dir/SuperLU/4.3/include
    cd $nalu_build_dir/packages/SuperLU_4.3
    make
    cp SRC/*.h $nalu_install_dir/SuperLU/4.3/include

If you are planning on sharing this install location with other team members, you will likely need to change a
permission on a particular file, ``chmod g+r  $nalu_install_dir/SuperLU/4.3/include/superlu_enum_consts.h``.

Libxml2 v2.9.2
~~~~~~~~~~~~~~

Libxml2 is found `here <http://www.xmlsoft.org/sources/>`__.

Prepare:

::

    cd $nalu_build_dir/packages
    curl -o libxml2-2.9.2.tar.gz http://www.xmlsoft.org/sources/libxml2-2.9.2.tar.gz
    tar -xvf libxml2-2.9.2.tar.gz

Build:

::

    cd $nalu_build_dir/packages/libxml2-2.9.2
    CC=mpicc CXX=mpicxx ./configure -without-python --prefix=$nalu_install_dir/libxml2/2.9.2
    make
    make install

Boost v1.68.0
~~~~~~~~~~~~~

Boost is found `here <http://www.boost.org>`__.

Prepare:

::

    cd $nalu_build_dir/packages
    curl -o boost_1_68_0.tar.gz http://iweb.dl.sourceforge.net/project/boost/boost/1.68.0/boost_1_68_0.tar.gz
    tar -zxvf boost_1_68_0.tar.gz

Build:

::

    cd $nalu_build_dir/packages/boost_1_68_0
    ./bootstrap.sh --prefix=$nalu_install_dir/boost/1.68.0 --with-libraries=signals,regex,filesystem,system,mpi,serialization,thread,program_options,exception

You may or may not need to edit ``project-config.jam`` and add a 'using mpi', e.g,

using mpi: /path/to/mpi/openmpi/bin/mpicc

::

    ./b2 -j 4 2>&1 | tee boost_build_one
    ./b2 -j 4 install 2>&1 | tee boost_build_intall

YAML-CPP 0.6.2
~~~~~~~~

YAML is provided `here <https://github.com/jbeder/yaml-cpp>`__. Versions of Nalu before v1.1.0 used earlier versions of YAML-CPP. For brevity only the 
latest build instructions are discussed and the history of the Nalu git repo can be used to find older installation instructions if required. YAML-CPP 
has introduced several fixes since v0.5.3 in the master branch, so it is recommended to use the 0.6.2 informal release.

Prepare:

::

    cd $nalu_build_dir/packages
    cd yaml-cpp 
    git checkout yaml-cpp-0.6.2

Build:

::

    cd $nalu_build_dir/packages/yaml-cpp
    mkdir build
    cd build
    cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CXX_FLAGS=-std=c++11 -DCMAKE_CC_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX=$nalu_install_dir/yaml/0.6.2 ..
    make
    make install


Zlib v1.2.11
~~~~~~~~~~~

Zlib is provided `here <http://www.zlib.net>`__.

Prepare:

::

    cd $nalu_build_dir/packages
    curl -o zlib-1.2.11.tar.gz http://zlib.net/zlib-1.2.11.tar.gz
    tar -zxvf zlib-1.2.11.tar.gz

Build:

::

    cd $nalu_build_dir/packages/zlib-1.2.11
    CC=gcc CXX=g++ CFLAGS=-O3 CXXFLAGS=-O3 ./configure --prefix=$nalu_install_dir/zlib/1.2.11
    make
    make install

HDF5 v1.10.4
~~~~~~~~~~~~

HDF5 1.10.4 is provided `here <http://www.hdfgroup.org/downloads/index.html>`__.

Prepare:

::

    cd $nalu_build_dir/packages/
    tar -zxvf hdf5-1.10.4.tar.gz

Build:

::

    cd $nalu_build_dir/packages/hdf5-1.10.4
    ./configure CC=mpicc FC=mpif90 CXX=mpicxx CXXFLAGS="-fPIC -O3" CFLAGS="-fPIC -O3" FCFLAGS="-fPIC -O3" --enable-parallel --with-zlib=$nalu_install_dir/zlib/1.2.11 --prefix=$nalu_install_dir/hdf5/1.10.4
    make
    make install
    make check

NetCDF v4.6.1 and Parallel NetCDF v1.10.1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to support all aspects of Nalu's parallel models, NetCDF and Parallel NetCFD must be consistent.

Parallel NetCDF v1.10.0
**********************

Parallel NetCDF is provided on the `Argon Trac Page <https://trac.mcs.anl.gov/projects/parallel-netcdf/wiki/Download>`__. Newer versions
can be found managed by `Northwestern <http://cucis.ece.northwestern.edu/projects/PnetCDF/download.html>`__.

Prepare:

::

    cd $nalu_build_dir/packages/
    tar -zxvf parallel-netcdf-1.10.0.tar.gz

Build:

::

    cd parallel-netcdf-1.10.0
    ./configure --prefix=$nalu_install_dir/pnetcdf/1.10.0 CC=mpicc FC=mpif90 CXX=mpicxx CFLAGS="-I$nalu_install_dir/pnetcdf/1.10.0/include -O3" LDFLAGS=-L$nalu_install_dir/pnetcdf/1.10.0/lib --disable-fortran
    make
    make install

Note that we have created an install directory that might look like ``$nalu_build_dir/install``.

NetCDF v4.6.1
***************

NetCDF is provided `here <https://github.com/Unidata/netcdf-c/releases>`__.

Prepare:

::

    cd $nalu_build_dir/packages/
    curl -o netcdf-c-4.6.1.tar.gz https://codeload.github.com/Unidata/netcdf-c/tar.gz/v4.6.1
    tar -zxvf netcdf-c-4.6.1.tar.gz 

Build:

::

    cd netcdf-c-4.6.1/
    ./configure --prefix=$nalu_install_dir/netcdf/4.6.1 CC=mpicc FC=mpif90 CXX=mpicxx CFLAGS="-I$nalu_install_dir/hdf5/1.10.4/include -I$nalu_install_dir/pnetcdf/1.10.0/include -O3" CPPFLAGS=${CFLAGS} LDFLAGS="-L$nalu_install_dir/hdf5/1.10.4/lib -L$nalu_install_dir/pnetcdf/1.10.0/lib -Wl,--rpath=$nalu_install_dir/hdf5/1.10.4/lib" --enable-pnetcdf --enable-parallel-tests --enable-netcdf-4 --disable-shared --disable-fsync --disable-cdmremote --disable-dap --disable-doxygen --disable-v2
    make -j 4 
    make check
    make install


Trilinos
~~~~~~~~

Trilinos is managed by the `Trilinos <http://www.trilinos.org>`__ project and can be found on Github.
The Nalu code base follows ``develop`` branch. 

Prepare:

::

    cd $nalu_build_dir/packages/
    git clone https://github.com/trilinos/Trilinos.git
    cd $nalu_build_dir/packages/Trilinos
    git checkout develop
    mkdir build

HYPRE
~~~~~

Nalu can use HYPRE solvers and preconditioners, especially for Pressure Poisson
solves. However, this dependency is optional and is not enabled by default.
Users wishing to use HYPRE solver and preconditioner combination must compile
HYPRE library and link to it when building Nalu.

.. code-block:: bash

   # 1. Clone hypre sources
   https://github.com/LLNL/hypre.git
   cd hypre/src

   # 2. Configure HYPRE package and pass installation directory
   ./configure --prefix=$nalu_install_dir --without-superlu --without-openmp --enable-bigint

   # 3. Compile and install
   make && make install

.. note::

   #. Make sure that ``--enable-bigint`` option is turned on if you intend to
      run linear systems with :math:`> 2` billion rows. Otherwise, ``nalu``
      executable will throw an error at runtime for large problems.

   #. Users must pass ``-DENABLE_HYPRE`` option to CMake during Nalu
      configuration phase. Optionally, the variable `-DHYPRE_DIR`` can be used
      to pass the path of HYPRE install location to CMake.

Build
*****

Place into the build directory, one of the ``do-configTrilinos_*`` files, that can be obtained from the Nalu repo.

``do-configTrilinos_*`` will be used to run cmake to build trilinos correctly for Nalu. Note that there are two files: one 
for 'release' and the other 'debug'. The files can be found on the Nalu GitHub site or copied from ``$nalu_build_dir/packages/Nalu/build``, 
which is created in the Nalu build step documented below. For example:

Pull latest version of ``do-configTrilinos_*`` from Nalu's GitHub site:

::

    curl -o $nalu_build_dir/packages/Trilinos/build/do-configTrilinos_release https://raw.githubusercontent.com/NaluCFD/Nalu/master/build/do-configTrilinos_release

Or if you create the Nalu directory as directed below, simply copy one of the ``do-configTrilinos_*`` files from local copy of Nalu's git repository:

::

    cp $nalu_build_dir/packages/Nalu/build/do-configTrilinos_release $nalu_build_dir/packages/Trilinos/build

Now edit ``do-configTrilinos_release`` to modify the paths so they point to the proper TPL ``$nalu_install_dir``.

::

    cd $nalu_build_dir/packages/Trilinos/build
    chmod +x do-configTrilinos_release

Make sure all other paths to netcdf, hdf5, etc., are correct.

::

    ./do-configTrilinos_release
    make
    make install


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

Nalu is provided `here <https://github.com/NaluCFD/Nalu>`__. One may either build the released Nalu version 1.2.0 which matches 
with Trilinos version 12.12.1, or the master branch of Nalu which matches with the master branch or develop branch of Trilinos. If 
it is necessary to build an older version of Nalu, refer to the history of the Nalu git repo for instructions on doing so.

Prepare:

::

    git clone https://github.com/NaluCFD/Nalu.git


Build
*****

In ``Nalu/build``, you will find the `do-configNalu <https://github.com/NaluCFD/Nalu/blob/master/build/do-configNalu_release>`__ script. 
Copy the ``do-configNalu_release`` or ``debug`` file to a new, non-tracked file:

::

    cp do-configNalu_release do-configNaluNonTracked

Edit the paths at the top of the files by defining the ``nalu_install_dir`` variable. Within ``Nalu/build``, execute the following commands:

::

    ./do-configNaluNonTracked
    make 

This process will create ``naluX`` within the ``Nalu/build`` location. You may also build a debug executable by modifying the Nalu 
config file to use "Debug". In this case, a ``naluXd`` executable is created.


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


