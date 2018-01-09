#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================

# Standard regression test
function(add_test_r testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r)

# Standard performance test
function(add_test_p testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 2500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "performance")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_p)

# Regression test with single restart
function(add_test_r_rst testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE} && mpiexec -np ${np} ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_rst.i -o ${testname}_rst.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname}_rst ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_rst.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r_rst)

# Verification test with three resolutions
function(add_test_v3 testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R0.i -o ${testname}_R0.log && mpiexec -np ${np} ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R1.i -o ${testname}_R1.log && mpiexec -np ${np} ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R2.i -o ${testname}_R2.log && python ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/norms.py")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "verification")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_v3)

# Verification test with two resolutions
function(add_test_v2 testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R0.i -o ${testname}_R0.log && mpiexec -np ${np} ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R1.i -o ${testname}_R1.log && python ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/norms.py")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "verification")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_v2)

# Regression test that runs with different numbers of processes
function(add_test_r_np testname np)
    add_test(${testname}Np${np} sh -c "mpiexec -np ${np} ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}Np${np}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname}Np${np} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}Np${np}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname}Np${np} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r_np)

# Standard unit test
function(add_test_u testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} ${CMAKE_BINARY_DIR}/${utest_ex_name}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1000 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "unit")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_u)

# Regression test with catalyst capability
function(add_test_r_cat testname np ncat)
    if(ENABLE_PARAVIEW_CATALYST)
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i.in)
        add_test(${testname} sh -c "mpiexec -np ${np} ${CMAKE_BINARY_DIR}/${nalu_ex_catalyst_name} -i ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}/${testname}_catalyst.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail_catalyst.sh ${testname} ${ncat}")
        set_tests_properties(${testname} PROPERTIES TIMEOUT 1000 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
        file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
        set(CATALYST_FILE_INPUT_DECK_COMMAND "catalyst_file_name: catalyst.txt")       
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i.in
                       ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}/${testname}_catalyst.i @ONLY)
        file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/catalyst.txt
             DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
      endif()
    else()
      add_test_r(${testname} ${np})
    endif()
endfunction(add_test_r_cat)

#=============================================================================
# Regression tests
#=============================================================================
add_test_r_cat(ablNeutralEdge 8 11)
add_test_r(ablStableElem 4)
add_test_r_rst(ablUnstableEdge 4)
add_test_r(ablUnstableEdge_ra 4)
add_test_r(actuatorLine 8)
add_test_r(concentricRad 4)
add_test_r(cvfemHC 8)
add_test_r(dgMMS 6)
add_test_r(dgNonConformal 4)
add_test_r(dgNonConformal3dFluids 4)
add_test_r(dgNonConformal3dFluidsHexTet 4)
add_test_r(dgNonConformal3dFluidsP1P2 8)
add_test_r(dgNonConformalEdge 4)
add_test_r(dgNonConformalEdgeCylinder 8)
add_test_r(dgNonConformalElemCylinder 8)
add_test_r(dgNonConformalFluids 4)
add_test_r(dgNonConformalFluidsEdge 4)
add_test_r_rst(dgNonConformalThreeBlade 4)
add_test_r(ductElemWedge 2)
add_test_r(ductWedge 2)
add_test_r(edgeHybridFluids 8)
add_test_r(edgePipeCHT 4)
add_test_r(ekmanSpiral 4)
add_test_r(ekmanSpiralConsolidated 4)
add_test_r(elemBackStepLRSST 4)
add_test_r(elemClosedDomain 2)
add_test_r(elemHybridFluids 8)
add_test_r(elemHybridFluidsShift 8)
add_test_r(elemPipeCHT 4)
add_test_r(femHC 2)
add_test_r(femHCGL 2)
add_test_r(fluidsPmrChtPeriodic 8)
add_test_r(heatedBackStep 4)
add_test_r_rst(heatedWaterChannelEdge 4)
add_test_r(heatedWaterChannelElem 4)
add_test_r_rst(heliumPlume 8)
add_test_r(hoHelium 8)
add_test_r(hoVortex 2)
add_test_r(inputFireEdgeUpwind 4)
add_test_r(inputFireElem 4)
add_test_r(kovasznay_P5 1)
add_test_r(milestoneRun 4)
add_test_r(milestoneRunConsolidated 4)
add_test_r_cat(mixedTetPipe 8 7)
add_test_r(movingCylinder 4)
add_test_r(nonConformalWithPeriodic 2)
add_test_r(nonConformalWithPeriodicConsolidated 2)
add_test_r(nonIsoEdgeOpenJet 4)
add_test_r(nonIsoElemOpenJet 4)
add_test_r(nonIsoElemOpenJetConsolidated 4)
add_test_r(nonIsoNonUniformEdgeOpenJet 4)
add_test_r(nonIsoNonUniformElemOpenJet 4)
add_test_r(overset 6)
add_test_r(oversetFluids 6)
add_test_r(oversetFluidsEdge 6)
add_test_r_np(periodic3dElem 1)
add_test_r_np(periodic3dElem 4)
add_test_r_np(periodic3dElem 8)
add_test_r_np(periodic3dEdge 1)
add_test_r_np(periodic3dEdge 4)
add_test_r_np(periodic3dEdge 8)
add_test_r(quad9HC 2)
add_test_r_cat(steadyTaylorVortex 4 6)
add_test_r(variableDensNonIso 2)
add_test_r(variableDensNonUniform 2)
add_test_r(variableDensNonUniform_P5 8)
if(ENABLE_OPENFAST)
   add_test_r(nrel5MWactuatorLine 4)
   add_subdirectory(test_files/nrel5MWactuatorLine)
endif(ENABLE_OPENFAST)
if(ENABLE_TIOGA)
  add_test_r(oversetSphereTIOGA 8)
endif(ENABLE_TIOGA)
if(ENABLE_HYPRE)
  add_test_r(dgncThreeBladeHypre 2)
endif(ENABLE_HYPRE)

#=============================================================================
# Convergence tests
#=============================================================================
add_test_v2(BoussinesqNonIso 8)
add_test_v3(cvfemHexHC_P3 8)
add_test_v3(hoVortex_P2 8)
add_test_v3(steadyTaylorVortex_P4 8)

#=============================================================================
# Unit tests
#=============================================================================
add_test_u(unitTest1 1)
add_test_u(unitTest2 2)

#=============================================================================
# Performance tests
#=============================================================================
add_test_p(oversetHybrid 8)
add_test_p(uqSlidingMeshDG 8)
add_test_p(waleElemXflowMixFrac3.5m 8)
