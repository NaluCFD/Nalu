# Testing with CTest

Testing is enabled through CTest and CDash and the results are posted
[here](http://my.cdash.org/index.php?project=Nalu).

## Prerequisites

- All third party libraries have been built (including Trilinos).

- The do-configNalu scripts point to the right directories. This
  should be changed at some point. Maybe SPACK can solve this.
    
## Running
	  
The CTest script can be run using the following command: 
```{bash}
cd Nalu/build
ctest -DNIGHTLY_DIR=$nightlyDirectory -S ../cmake/$ctest_script
```
where `$nightlyDirectory` is the directory where the tests are run.
