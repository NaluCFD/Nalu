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
cd build
# ./do-configTrilinos
# make -j 4
make
make install >> $baseTestDirectory/nightlyBuildAndTest/TrilinosBuild.txt

# get to Nalu
cd $baseTestDirectory/nightlyBuildAndTest/Nalu

# pull latest Nalu
git pull

# get to build
cd build

# config
./do-ConfigNaluNonTracked

# build it... send contents to a file
make -j 4 >> $baseTestDirectory/nightlyBuildAndTest/NaluBuild.txt

# get to NaluRtest; do not hold an extra clone - switch to nightly
cd  $baseTestDirectory/NaluRtest
git checkout nightly

# pull
git pull

# remove old test vouchers
if [ -d "$baseTestDirectory/runNaluRtest" ]; then
    rm -rf $baseTestDirectory/runNaluRtest/nightly/*/PASS
fi

# run it... send contents to a file
./run_tests.sh >> $baseTestDirectory/nightlyBuildAndTest/NaluRtest.txt

# checkout master
git checkout master

# get back to base
cd $baseTestDirectory/nightlyBuildAndTest

# mail contents
osascript "mailRtest.scpt"
