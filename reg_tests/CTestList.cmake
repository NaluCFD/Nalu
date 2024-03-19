#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================

# Standard regression test
function(add_test_r testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r)

# Standard performance test
function(add_test_p testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 2500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "performance")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_p)

# Regression test with single restart
function(add_test_r_rst testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE} && mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_rst.i -o ${testname}_rst.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname}_rst ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_rst.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r_rst)

# Standard regression test with input data
function(add_test_r_t testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r_t)

# Verification test with three resolutions
function(add_test_v3 testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R0.i -o ${testname}_R0.log && mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R1.i -o ${testname}_R1.log && mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R2.i -o ${testname}_R2.log && python ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/norms.py")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "verification")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_v3)

# Verification test with two resolutions
function(add_test_v2 testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R0.i -o ${testname}_R0.log && mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R1.i -o ${testname}_R1.log && python ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/norms.py")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "verification")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_v2)

# Regression test that runs with different numbers of processes
function(add_test_r_np testname np)
    add_test(${testname}Np${np} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}Np${np}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname}Np${np} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}Np${np}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname}Np${np} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r_np)

# Standard unit test
function(add_test_u testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${utest_ex_name}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1000 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "unit")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_u)

# Regression test with catalyst capability
function(add_test_r_cat testname np ncat cat_comm cat_file)
    if(ENABLE_PARAVIEW_CATALYST)
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i.in)
        add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_catalyst_name} -i ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}/${testname}_catalyst.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail_catalyst.sh ${testname} ${ncat}")
        set_tests_properties(${testname} PROPERTIES TIMEOUT 1000 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
        file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
        set(CATALYST_FILE_INPUT_DECK_COMMAND "${cat_comm} ${cat_file}")
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i.in
                       ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}/${testname}_catalyst.i @ONLY)
        file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${cat_file}
             DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
      endif()
    else()
      add_test_r(${testname} ${np})
    endif()
endfunction(add_test_r_cat)

function(add_test_l testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/laboratory/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/laboratory/${testname}/${testname}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 2500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/laboratory/${testname}" LABELS "laboratory")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/laboratory/${testname})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/laboratory/${testname}/mesh)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_files/laboratory/${testname}/mesh/${testname}.exo DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test_files/laboratory/${testname}/mesh)
endfunction(add_test_l)

#=============================================================================
# Regression tests
#=============================================================================
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
add_test_r_t(ductWedge 2)
add_test_r(edgeHybridFluids 8)
add_test_r(edgePipeCHT 4)
add_test_r_rst(elemBackStepLRSST 4)
add_test_r(elemClosedDomain 2)
add_test_r(elemHybridFluids 8)
add_test_r(elemHybridFluidsShift 8)
add_test_r(elemPipeCHT 4)
add_test_r(femHC 2)
add_test_r(femHCGL 2)
add_test_r(fluidsPmrChtPeriodic 8)
add_test_r_rst(fluidsPmrChtPeriodic 8)
add_test_r(heatedBackStep 4)
add_test_r(1x2x10Tet10 4)
add_test_r(ductTet10 4)
add_test_r_rst(heatedWaterChannelEdge 4)
add_test_r(heatedWaterChannelElem 4)
add_test_r_rst(heliumPlume 8)
add_test_r(hoHelium 8)
add_test_r(hoVortex 2)
add_test_r(inputFireEdgeUpwind 4)
add_test_r(inputFireElem 4)
add_test_r(milestoneRun 4)
add_test_r(milestoneRunConsolidated 4)
add_test_r_cat(mixedTetPipe 8 7 "paraview_script_name:" "pipe_tet_catalyst.py")
add_test_r(movingCylinder 4)
add_test_r(nonConformalWithPeriodic 2)
add_test_r(nonConformalWithPeriodicConsolidated 2)
add_test_r(nonIsoEdgeOpenJet 4)
add_test_r(nonIsoElemOpenJetConsolidated 4)
add_test_r(nonIsoNonUniformEdgeOpenJet 4)
add_test_r(nonIsoNonUniformElemOpenJet 4)
add_test_r(nonUniformElemOpenJet 4)
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
add_test_r(sphereDrop 4)
add_test_r_cat(steadyTaylorVortex 4 6 "catalyst_file_name:" "catalyst.txt")
add_test_r(variableDensNonIso 2)
add_test_r(variableDensNonUniform 2)
add_test_r(femPassiveScalar 4)
add_test_r(femFluidsVortex 4)
add_test_r_rst(kEpsPipe 8)
add_test_r(vofBuoy 4)
add_test_r(vofSlosh 4)
add_test_r(vofWaveGenerator 4)
add_test_r(dgNonConformalShockTube 2)

#=============================================================================
# Convergence tests
#=============================================================================
add_test_v2(BoussinesqNonIso 8)

#=============================================================================
# Unit tests
#=============================================================================
add_test_u(unitTest1 1)
add_test_u(unitTest2 2)

#=============================================================================
# Performance tests
#=============================================================================
if(ENABLE_PERFORMANCE_TESTS)
  add_test_p(oversetHybrid 8)
  add_test_p(uqSlidingMeshDG 8)
  add_test_p(waleElemXflowMixFrac3.5m 8)
endif(ENABLE_PERFORMANCE_TESTS)

#=============================================================================
# Laboratory tests
#=============================================================================
if(ENABLE_LABORATORY_TESTS)
  add_test_l(2d_quad4_channel 2)
  add_test_l(2d_quad9_couette 4)
  add_test_l(3d_hex8_open_jet 8)
  add_test_l(3d_tet4_pipe 8)
  add_test_l(3d_hex8_dam_break 8)
  add_test_l(2d_quad9_helium 8)
  add_test_l(2d_tri3_quad4_street 8)
  add_test_l(1d_quad4_adv_diff 1)
  add_test_l(3d_tet4_taylor_green_0p2 8)
  add_test_l(3d_hex8_turb_channel 4)
endif(ENABLE_LABORATORY_TESTS)
