#!/bin/bash

# assign base directory
baseTestDirectory=/Users/naluIt/gitHubWork

# get date
currentDate="$(date +'%d/%m/%Y')"
echo $currentDate

# get to nightly test directory
cd $baseTestDirectory/nightlyBuildAndTest

# assign files to current date
echo $currentDate > NaluBuild.txt
echo $currentDate > NaluRtest.txt
echo $currentDate > TrilinosBuild.txt

# get to Trilinos and pull
cd $baseTestDirectory/packages/publicTrilinos
#git pull

# now build
echo ....Commencing the Trilinos Build.....
cd build
./do-configTrilinos
make -j 4
make install >> $baseTestDirectory/nightlyBuildAndTest/TrilinosBuild.txt

# get to Nalu
cd $baseTestDirectory/nightlyBuildAndTest/Nalu

# pull latest Nalu
git pull

# get to build
cd build

# remove old executable
if [ -f "naluX" ]; then
        echo Removing previous naluX executable
	rm naluX
fi

# config
./do-ConfigNaluNonTracked >> $baseTestDirectory/nightlyBuildAndTest/NaluBuild.txt

# build it... send contents to a file
echo ....Commencing the Nalu Build....
make -j 4 >> $baseTestDirectory/nightlyBuildAndTest/NaluBuild.txt

# check for executable
if [ -f "naluX" ]; then
        echo Nalu build successful... ready to test..
	# get to NaluRtest; do not hold an extra clone - switch to nightly
	cd  $baseTestDirectory/NaluRtest
	git checkout nightly

	# pull
	git pull

	# remove old test vouchers
	if [ -d "$baseTestDirectory/runNaluRtest" ]; then
    		rm -rf $baseTestDirectory/runNaluRtest/nightly/*/PASS
    		echo ....Removing PASS status under runNaluRtest....
	fi

	# run it... send contents to a file
	echo ....Commencing NaluRtest....
	./run_tests.sh >> $baseTestDirectory/nightlyBuildAndTest/NaluRtest.txt

	# checkout master
	git checkout master

else
	echo PROCESS FAILED >> $baseTestDirectory/nightlyBuildAndTest/NaluBuild.txt
        echo PROCESS FAILED >> $baseTestDirectory/nightlyBuildAndTest/NaluRtest.txt	
fi

# get back to base
cd $baseTestDirectory/nightlyBuildAndTest

# mail contents
echo ....Mailing the results voucher....
osascript "mailRtest.scpt"
