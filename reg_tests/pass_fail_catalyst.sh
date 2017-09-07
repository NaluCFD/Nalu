#!/bin/bash

determine_pass_fail() {

    for file in catalyst_test_image_output/CatalystTestImage*.png; do 
        if [ -s ${file} ]; then
            let foundImageCount=foundImageCount+1
        fi
    done

    if [ ${foundImageCount} -ne ${numImages} ]; then 
        return 1
    else
        return 0
    fi

}

# $1 is the test name
main() {
  testName=$1
  numImages=$2
  foundImageCount=0
  determine_pass_fail ${numImages}
  passStatus="$?"
  performanceTime=`grep "STKPERF: Total Time" ${testName}.log  | awk '{print $4}'`
  if [ ${passStatus} -eq 0 ]; then
      echo -e "..${testName} Catalyst........... PASSED":" " ${performanceTime} " s"
      exit 0
  else
      echo -e -n "..${testName} Catalyst........... FAILED":" " ${performanceTime} " s"
      echo -e ", expected " ${numImages} " images, and found " ${foundImageCount} " images"
      exit 1
  fi
}

# ./pass_fail_catalyst.sh testName numImages
main "$@"
