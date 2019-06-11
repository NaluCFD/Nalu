#!/bin/bash

# make sure permissions are set
kinit -k -t /ascldap/users/spdomin/private/spdomin.keytab spdomin

# assign base directory
moduleDirectory=/ascldap/users/spdomin/gitHubWork/modules
baseTrilinosDirectory=/home/spdomin/gitHubWork/scratch_build/packages/Trilinos

# make sure our proxy's are set
cd $moduleDirectory
source gitHubProxy.txt
source nalu_module_gcc7.2.0_mpi1.10.2.sh

# leave the pool
leavePool

# get date
currentDate="$(date +'%m_%d_%Y')"
echo $currentDate

TrilinosBuildFileName="TrilinosStableUpdate_$currentDate.txt"

# get to Trilinos and pull on stable
cd $baseTrilinosDirectory

echo "$currentDate Trilinos Stable Update" >> $baseTrilinosDirectory/$TrilinosBuildFileName

echo "....Checkout stable Trilinos..."  >> $baseTrilinosDirectory/$TrilinosBuildFileName

git checkout stable 
git pull origin develop 

# delete old stuff
cp build_stable_release_7.2.0/do-configTrilinos_stable_release_7.2.0 ./
echo "....Deleting Trilinos Stable Release build directory...." >> $baseTrilinosDirectory/$TrilinosBuildFileName
rm -rf build_stable_release_7.2.0

# re-construct build
mkdir build_stable_release_7.2.0
mv do-configTrilinos_stable_release_7.2.0 build_stable_release_7.2.0

# now build
echo "....Commencing the Trilinos Stable Release Build...." >> $baseTrilinosDirectory/$TrilinosBuildFileName
cd build_stable_release_7.2.0

echo "....do-configTrilinos Stable Release" >> $baseTrilinosDirectory/$TrilinosBuildFileName
./do-configTrilinos_stable_release_7.2.0

echo "....make on Trilinos Stable Release"  >> $baseTrilinosDirectory/$TrilinosBuildFileName
make -j 8
make install
echo "....install on Trilinos Stable Release Complete" >> $baseTrilinosDirectory/$TrilinosBuildFileName

# now debug....
cd $baseTrilinosDirectory

# delete old stuff
cp build_stable_debug_7.2.0/do-configTrilinos_stable_debug_7.2.0 ./
echo "....Deleting Trilinos Stable Debug build directory...." >> $baseTrilinosDirectory/$TrilinosBuildFileName
rm -rf build_stable_debug_7.2.0

# re-construct build
mkdir build_stable_debug_7.2.0
mv do-configTrilinos_stable_debug_7.2.0 build_stable_debug_7.2.0

# now build
echo "....Commencing the Trilinos Stable Debug Build...." >> $baseTrilinosDirectory/$TrilinosBuildFileName
cd build_stable_debug_7.2.0

echo "....do-configTrilinos Stable Debug" >> $baseTrilinosDirectory/$TrilinosBuildFileName
./do-configTrilinos_stable_debug_7.2.0

echo "....make on Trilinos Stable Debug" >> $baseTrilinosDirectory/$TrilinosBuildFileName
make -j 8
make install
echo "....install on Trilinos Stable Debug Complete" >> $baseTrilinosDirectory/$TrilinosBuildFileName

# get back to base Trilinos
cd $baseTrilinosDirectory 

# mail contents; just the tests for now...
echo "....Mailing the results voucher...." >> $baseTrilinosDirectory/$TrilinosBuildFileName 
cat $TrilinosBuildFileName | mail -s "TrilinosStableUpdate; Linux: $currentDate" spdomin@sandia.gov 


