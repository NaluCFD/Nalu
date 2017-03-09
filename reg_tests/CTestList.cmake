
function(add_test_r testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} ${CMAKE_BINARY_DIR}/naluX -i ${testname}.i -o ${testname}.log && ${CMAKE_BINARY_DIR}/pass_fail.sh ${testname} ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 500 PROCESSORS ${np} WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/${testname}")
endfunction(add_test_r)

function(add_test_u testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} ${CMAKE_BINARY_DIR}/unittestX")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 500 PROCESSORS ${np} WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/unitTests")
endfunction(add_test_u)

function(add_test_r_rst testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} ${CMAKE_BINARY_DIR}/naluX -i ${testname}.i -o ${testname}.log && ${CMAKE_BINARY_DIR}/pass_fail.sh ${testname} ${TOLERANCE} && mpiexec -np ${np} ${CMAKE_BINARY_DIR}/naluX -i ${testname}_rst.i -o ${testname}_rst.log && ${CMAKE_BINARY_DIR}/pass_fail.sh ${testname}_rst ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 500 PROCESSORS ${np} WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/${testname}")
endfunction(add_test_r_rst)

function(add_test_r_np testname np)
    add_test(${testname}Np${np} sh -c "mpiexec -np ${np} ${CMAKE_BINARY_DIR}/naluX -i ${testname}.i -o ${testname}Np${np}.log && ${CMAKE_BINARY_DIR}/pass_fail.sh ${testname}Np${np} ${TOLERANCE}")
    set_tests_properties(${testname}Np${np} PROPERTIES TIMEOUT 500 PROCESSORS ${np} WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/${testname}")
endfunction(add_test_r_np)

function(add_test_p testname np timeo)
    add_test(${testname} sh -c "mpiexec -np ${np} ${CMAKE_BINARY_DIR}/naluX -i ${testname}.i -o ${testname}.log && ${CMAKE_BINARY_DIR}/pass_fail.sh ${testname} ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT ${timeo} PROCESSORS ${np} WORKING_DIRECTORY "${PERF_TEST_RESULT_DIRECTORY}/${testname}")
endfunction(add_test_p)

#=============================================================================
# Regression tests
#=============================================================================

set(this_test periodic3dElem)
add_test_r_np(${this_test} 1)
add_test_r_np(${this_test} 4)
add_test_r_np(${this_test} 8)
set(this_test periodic3dEdge)
add_test_r_np(${this_test} 1)
add_test_r_np(${this_test} 4)
add_test_r_np(${this_test} 8)

add_test_r(quad9HC 2)
add_test_r(steadyTaylorVortex 4)
add_test_r(hoVortex 2)
add_test_r(hoHelium 8)
add_test_r(dgNonConformal 4)
add_test_r(dgNonConformalEdge 4)
add_test_r(dgNonConformalFluids 4)
add_test_r(dgNonConformalFluidsEdge 4)
add_test_r(dgNonConformal3dFluids 4)
add_test_r(dgNonConformal3dFluidsP1P2 8)
add_test_r(dgNonConformal3dFluidsHexTet 4)
add_test_r_rst(dgNonConformalThreeBlade 4)
add_test_r(overset 6)
add_test_r(oversetFluids 6)
add_test_r(oversetFluidsEdge 6)
add_test_r(concentricRad 4)
add_test_r(movingCylinder 4)
add_test_r(elemBackStepLRSST 4)
add_test_r(ductWedge 2)
add_test_r(ductElemWedge 2)
add_test_r(heatedBackStep 4)
add_test_r(edgePipeCHT 4)
add_test_r(elemPipeCHT 4)
#add_test_r_rst(heliumPlume 8) # Input file does not exist error
add_test_r(dgNonConformalEdgeCylinder 8)
add_test_r(dgNonConformalElemCylinder 8)
add_test_r(fluidsPmrChtPeriodic 8)
add_test_r(nonIsoElemOpenJet 4)
add_test_r(nonIsoEdgeOpenJet 4)
#add_test_r(hdf5VarZChi 4) # This one was already disabled
add_test_r(elemHybridFluids 8)
add_test_r(elemHybridFluidsShift 8)
add_test_r(edgeHybridFluids 8)
add_test_r(elemClosedDomain 2)
add_test_r(mixedTetPipe 8)
add_test_r(inputFireEdgeUpwind 4)
add_test_r(inputFireElem 4)
add_test_r(nonIsoNonUniformElemOpenJet 4)
add_test_r(nonIsoNonUniformEdgeOpenJet 4)
add_test_r(milestoneRun 4)
#add_test_r_rst(heatedWaterChannel 4) # Maybe need to rewrite this one
#add_test_r(variableDensMMS 2) # Input file does not exist error
add_test_r(actuatorLine 8)
add_test_r(femHC 2)
add_test_r(cvfemHC 8)
add_test_r(ablUnstableEdge 4)
add_test_r(ablStableElem 4)
add_test_r(ekmanSpiral 4)
add_test_r(dgMMS 6)

#=============================================================================
# Unit tests
#=============================================================================

add_test_u(unitTest1 1)
add_test_u(unitTest2 2)

#=============================================================================
# Performance tests
#=============================================================================

#add_test_p(waleElemXflowMixFrac3 8 600) # Gives BAD_COMMAND error
add_test_p(uqSlidingMeshDG 8 1000)
add_test_p(oversetHybrid 8 1000)
