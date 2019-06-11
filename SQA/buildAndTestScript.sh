#!/bin/bash

# master or develop?
testMaster='0'

# make sure permissions are set
kinit spdomin@dce.sandia.gov -k -t /ascldap/users/spdomin/private/spdomin.keytab

# source 
source /ascldap/users/spdomin/.profile

# assign base directory
moduleDirectory=/ascldap/users/spdomin/gitHubWork/modules
baseTestDirectory=/scratch/spdomin

# make sure our proxy's are set
cd $moduleDirectory
source gitHubProxy.txt
source nalu_module_gcc7.2.0_mpi1.10.2.sh

cmake --version

# leave the pool
leavePool

# get date
currentDate="$(date +'%m/%d/%Y')"
echo $currentDate

# get to nightly test directory
cd $baseTestDirectory/nightlyBuildAndTest

# assign files to current date
echo $currentDate NaluBuild > NaluBuild.txt
echo $currentDate NaluRtest > NaluRtest.txt
echo $currentDate TrilinosBuild > TrilinosBuild.txt

# get to Trilinos and pull on master
cd Trilinos
if [ "$testMaster" == '1' ]; then
    echo ....Checkout master Trilinos... 
    git checkout master 
    git pull origin master
else
    echo ....Checkout develop Trilinos... 
    git checkout develop
    git pull origin develop 
fi

# get current SHA1
currentTrilinosSHA1="$(git log -1 --format="%H")"

# delete old stuff
cp build_nightly_release/do-configTrilinos_nightly_release ./
#echo ....Deleting Trilinos build directory....
rm -rf build_nightly_release
rm -rf /home/spdomin/gitHubWork/scratch_build/install/gcc7.2.0/Trilinos_nightly_release

# re-construct build
mkdir build_nightly_release
mv do-configTrilinos_nightly_release build_nightly_release

# now build
echo ....Commencing the Trilinos Build....
cd build_nightly_release

#echo ....do-configTrilinos_nightly_release
./do-configTrilinos_nightly_release

#echo ....make on Trilinos
make -j 16 >> $baseTestDirectory/nightlyBuildAndTest/TrilinosBuild.txt
make install >> $baseTestDirectory/nightlyBuildAndTest/TrilinosBuild.txt

# get to Nalu
cd $baseTestDirectory/nightlyBuildAndTest/NaluNightly

# pull latest Nalu
git checkout master
echo ....pull Nalu on origin master
git pull origin master

# update submodules
echo ...updating submodules
git submodule update

# get current SHA1
currentNaluSHA1="$(git log -1 --format="%H")"

# get to build
cd build

# delete old stuff
./cleanIt.sh

# remove old executable
if [ -f "naluX" ]; then
        echo Removing previous naluX executable
	rm naluX
fi

# config
echo ....do-configNalu_nightly_release
./do-configNalu_nightly_release >> $baseTestDirectory/nightlyBuildAndTest/NaluBuild.txt

# build it... send contents to a file
echo ....Commencing the Nalu Build....
make -j 16 >> $baseTestDirectory/nightlyBuildAndTest/NaluBuild.txt

# check for executable
if [ -f "naluX" ]; then
        echo Nalu build successful... ready to test..
	# get to nightlyBuildAndTest/NaluRtest; switch to nightly
	cd  $baseTestDirectory/nightlyBuildAndTest/NaluNightly/build

	# run it... send contents to a file
	echo ....Commencing NaluRtest....
        echo $pwd
	ctest -j 16 --output-on-failure >> $baseTestDirectory/nightlyBuildAndTest/NaluRtest.txt
else
	echo PROCESS FAILED >> $baseTestDirectory/nightlyBuildAndTest/NaluBuild.txt
        echo PROCESS FAILED >> $baseTestDirectory/nightlyBuildAndTest/NaluRtest.txt	
fi

# get back to base
cd $baseTestDirectory/nightlyBuildAndTest

# provide SHA1 info
echo NaluCFD/Nalu SHA1: $currentNaluSHA1 >> $baseTestDirectory/nightlyBuildAndTest/NaluRtest.txt
if [ "$testMaster" == '1' ]; then
    echo Trilinos/master SHA1: $currentTrilinosSHA1 >> $baseTestDirectory/nightlyBuildAndTest/NaluRtest.txt
else
    echo Trilinos/develop SHA1: $currentTrilinosSHA1 >> $baseTestDirectory/nightlyBuildAndTest/NaluRtest.txt
fi


# mail contents; just the tests for now...
echo ....Mailing the results voucher....
cat NaluRtest.txt | mail -s "NaluCFD/Nalu Build/Test Linux: $currentDate" spdomin@sandia.gov 
cat NaluRtest.txt | mail -s "NaluCFD/Nalu Build/Test Linux: $currentDate" nalucfd@gmail.com