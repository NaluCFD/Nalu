Building Nalu Semi-Automatically Using Spack
============================================

Mac OS X or Linux
-----------------

The following describes how to build Nalu and its dependencies
mostly automatically on your Mac using 
`Spack <https://spack.readthedocs.io/en/latest>`__. 
This can also be used as a template to build Nalu on any 
Linux system with Spack.

Step 1
~~~~~~

This assumes you have a Homebrew installation of GCC installed already 
(we are using GCC 6.3.0). These instructions have been tested on OSX 10.11. MacOS 10.12 
will not build CMake with GCC anymore, so these instructions won't work 
in that case, but we have built Nalu using Spack on MacOS Sierra by
using Homebrew to install ``cmake`` and ``pkg-config`` and defining these 
as external packages in Spack (see 
`packages.yaml.mac_sierra <https://github.com/NaluCFD/NaluSpack/blob/master/spack_config/packages.yaml.mac_sierra>`__).

Step 2
~~~~~~

Checkout the official Spack repo from github (we will checkout into ``${HOME}``):

``cd ${HOME} && git clone https://github.com/LLNL/spack.git``

Step 3
~~~~~~

Add Spack shell support to your ``.profile`` by adding the lines:

::

    export SPACK_ROOT=${HOME}/spack
    . $SPACK_ROOT/share/spack/setup-env.sh

Step 4
~~~~~~

Copy the ``nalu`` and ``nalu-trilinos`` package files from the 
`NaluSpack <https://github.com/NaluCFD/NaluSpack>`__ repo to
your installation of Spack or run the 
`copy_config.sh <https://github.com/NaluCFD/NaluSpack/blob/master/spack_config/copy_config.sh>`__
script from the repo:

::

    cd ${HOME} && git clone https://github.com/NaluCFD/NaluSpack.git
    cp -R ${HOME}/NaluSpack/spack_config/nalu ${SPACK_ROOT}/var/spack/repos/builtin/packages/
    cp -R ${HOME}/NaluSpack/spack_config/nalu-trilinos ${SPACK_ROOT}/var/spack/repos/builtin/packages/

Step 5
~~~~~~

Try ``spack info nalu`` to see if Spack works. If it does, check the
compilers you have available by:

::

    machine:~ user$ spack compilers
    ==> Available compilers
    -- gcc ----------------------------------------------------------
    gcc@6.3.0  gcc@4.2.1

    -- clang --------------------------------------------------------
    clang@8.0.0-apple  clang@7.3.0-apple

Step 6
~~~~~~

Install Nalu with whatever version of GCC (6.3.0 for us) you have
installed from Homebrew and force the build to use CMake 3.6.1 because
newer versions currently don't build on OS X:

::

    spack install nalu %gcc@6.3.0 ^cmake@3.6.1

That should be it! Spack will automatically use the most up-to-date dependencies 
unless otherwise specified. For example to constrain Nalu to use some specific 
versions of dependencies you could issue the Spack install command:

::

    spack install nalu %gcc@6.3.0 ^openmpi@1.10.3 ^boost@1.55.0 ^hdf5@1.8.16 ^parallel-netcdf@1.6.1 ^netcdf@4.3.3.1 ^cmake@3.6.1


NREL's Peregrine Machine
------------------------

The following describes how to build Nalu and its dependencies
mostly automatically on NREL's Peregrine machine using Spack. This can also be
used as a template to help build Nalu on any Linux system with Spack.

Step 1
~~~~~~

Login to Peregrine, and checkout the ``https://github.com/NaluCFD/NaluSpack.git`` 
repo (we will be cloning into the ${HOME} directory):

::

   cd ${HOME} && git clone https://github.com/NaluCFD/NaluSpack.git

One first thing to note is that the login nodes and the compute nodes on Peregrine 
run different OSes. So programs will be organized in Spack according to the OS 
they were built on, i.e. a login node (rhel6) typically called the front-end or 
compute node (centos6) typically called the back-end. You can see this in the 
directory structure where the programs will be built which will be located 
in ``${SPACK_ROOT}/opt``. You should build on a compute node.

Step 2
~~~~~~

Checkout the official Spack repo from github:

``cd ${HOME} && git clone https://github.com/LLNL/spack.git``

Step 3
~~~~~~

Configure your environment in the recommended way. You should purge all 
modules and only load GCC 5.2.0 in your login script. In the example 
`.bash_profile <https://github.com/NaluCFD/NaluSpack/blob/master/spack_config/dot_bash_profile.sh>`__
in the repo we also load Python. If you have problems building with Spack on 
Peregrine, it is most likely your environment has deviated from this 
recommended one. Even when building with the Intel compiler in Spack, 
this is the recommended environment.

::

   {
   module purge
   module load gcc/5.2.0
   module load python/2.7.8
   } &> /dev/null

Also add Spack shell support to your ``.bash_profile`` as shown in the example 
`.bash_profile <https://github.com/NaluCFD/NaluSpack/blob/master/spack_config/dot_bash_profile.sh>`__
in the repo or the following lines:

::

   export SPACK_ROOT=${HOME}/spack
   . $SPACK_ROOT/share/spack/setup-env.sh

Step 4
~~~~~~

Configure Spack for Peregrine. This is done by copying the ``compilers.yaml``, 
``config.yaml``, ``packages.yaml``, ``openmpi``, ``nalu``, and ``nalu-trilinos`` 
files/directories from the NaluSpack repo into your local ``${SPACK_ROOT}`` directory. 
These provide local configurations we need for Peregrine that override Spack's 
default configuration and the custom package files to install Nalu and the custom 
Trilinos build for Nalu. You can do this using the
`copy_config.sh <https://github.com/NaluCFD/NaluSpack/blob/master/spack_config/dot_bash_profile.sh>`__
script provided or by doing it manually as such:

::

   cp config.yaml ${SPACK_ROOT}/etc/spack/
   cp packages.yaml ${SPACK_ROOT}/etc/spack/
   cp compilers.yaml ${SPACK_ROOT}/etc/spack/
   cp -R nalu ${SPACK_ROOT}/var/spack/repos/builtin/packages/
   cp -R nalu-trilinos ${SPACK_ROOT}/var/spack/repos/builtin/packages/
   cp -R openmpi ${SPACK_ROOT}/var/spack/repos/builtin/packages/

Step 5
~~~~~~

Log out and log back in or source your ``.bash_profile`` to get the Spack 
shell support loaded. Try ``spack info nalu`` to see if Spack works.

Step 6
~~~~~~

Install Nalu using a compute node either interactively 
(``qsub -V -I -l nodes=1:ppn=24,walltime=4:00:00 -A <allocation> -q short``) 
or with the example batch script  
`install_nalu_gcc.sh <https://github.com/NaluCFD/NaluSpack/blob/master/install_scripts/install_nalu_gcc.sh>`__
(``qsub install_nalu_gcc.sh``):

::

   spack install binutils %gcc
   . ${SPACK_ROOT}/share/spack/setup-env.sh
   spack load binutils
   spack install nalu %gcc ^openmpi+tm+mxm@1.10.3

That's it! Hopefully the ``spack install nalu %gcc ^openmpi+tm+mxm@1.10.3`` 
command installs the entire set of dependencies and you get a working build 
of Nalu on Peregrine...after about 2 hours of waiting for it to build.

To build with the Intel compiler, note the necessary commands in 
`install_nalu_intel.sh <https://github.com/NaluCFD/NaluSpack/blob/master/install_scripts/install_nalu_intel.sh>`__ 
batch script (note you will need to point ``${TMPDIR}`` to disk as it defaults to 
RAM and will cause problems when building Trilinos, and also enabling accelerated 
message protocols segfaults the MPI wrappers when using the Intel compiler so 
remove ``+mxm`` too; also the latest m4 doesn't build with Intel).

Then to load Nalu (and you will need Spack's openmpi for Nalu now) into your path you 
will need to ``spack load openmpi %compiler`` and ``spack load nalu %compiler``, using 
``%gcc`` or ``%intel`` to specify which to load.

Development Build
-----------------

When building Nalu with Spack, Spack will cache downloaded archive files such as 
``*.tar.gz`` files. However, by default Spack will also erase extracted or
checked out ('staged') source files after it has built a package successfully. 
Therefore if your build succeeds, Spack will have erased the Nalu source code 
it checked out from Github. To get a version of Nalu you can develop in, 
the idea is to checkout Nalu yourself outside of Spack and build this version 
using the dependencies Spack has built for you. To do so, checkout Nalu:

::

   git clone https://github.com/NaluCFD/Nalu.git

Next, create your own directory, or use the existing ``build`` directory in Nalu to 
run the CMake configuration. When running the CMake configuration, point Nalu to 
the dependencies by using ``spack location -i <package>``. For example in the 
``build`` directory run:

::

   cmake -DTrilinos_DIR:PATH=`spack location -i nalu-trilinos` \
         -DYAML_DIR:PATH=`spack location -i yaml-cpp` \
         -DCMAKE_BUILD_TYPE=RELEASE \
         ..
   make

This should allow you to have a build of Nalu in which you are able to continuosly modify 
the source code and rebuild.

One more thing to note is there is an option for Trilinos in which you can build with 
debug symbols. At the ``spack install`` command, to turn debug on (we use the CMake ``RelWithDebInfo`` 
type), specify it with ``^nalu-trilinos+debug@master``. Then in your development build 
of Nalu, you can specify the ``RelWithDebInfo`` build type in CMake.
