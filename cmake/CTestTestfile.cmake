# CMake Testfile configures all the tests that we will run
# These tests directories are prepared in prep_tests.cmake


#=============================================================================
#
# Nightly tests
#
#=============================================================================

#=============================================================================
# periodic3dElem test
#=============================================================================
add_test(periodic3dElemNp1 "periodic3dElemNp1.sh")
set_tests_properties(periodic3dElemNp1
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 1
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/periodic3dElem")

add_test(periodic3dElemNp4 "periodic3dElemNp4.sh")
set_tests_properties(periodic3dElemNp4
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/periodic3dElem")

add_test(periodic3dElemNp8 "periodic3dElemNp8.sh")
set_tests_properties(periodic3dElemNp8
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 8
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/periodic3dElem")


#=============================================================================
# periodic3dEdge test
#=============================================================================
add_test(periodic3dEdgeNp1 "periodic3dEdgeNp1.sh")
set_tests_properties(periodic3dEdgeNp1
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 1
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/periodic3dEdge")

add_test(periodic3dEdgeNp4 "periodic3dEdgeNp4.sh")
set_tests_properties(periodic3dEdgeNp4
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/periodic3dEdge")

add_test(periodic3dEdgeNp8 "periodic3dEdgeNp8.sh")
set_tests_properties(periodic3dEdgeNp8
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 8
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/periodic3dEdge")


#=============================================================================
# quad9HC test
#=============================================================================
add_test(quad9HC "quad9HC.sh")
set_tests_properties(quad9HC
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 2
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/quad9HC")


#=============================================================================
# steadyTaylorVortex test
#=============================================================================
add_test(steadyTaylorVortex "steadyTaylorVortex.sh")
set_tests_properties(steadyTaylorVortex
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/steadyTaylorVortex")


#=============================================================================
# hoVortex test
#=============================================================================
add_test(hoVortex "hoVortex.sh")
set_tests_properties(hoVortex
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 1000
              PROCESSORS 2
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/hoVortex")


#=============================================================================
# hoHelium test
#=============================================================================
add_test(hoHelium "hoHelium.sh")
set_tests_properties(hoHelium
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 600
              PROCESSORS 8
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/hoHelium")


#=============================================================================
# dgNonConformal test
#=============================================================================
add_test(dgNonConformal "dgNonConformal.sh")
set_tests_properties(dgNonConformal
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformal")


#=============================================================================
# dgNonConformalEdge test
#=============================================================================
add_test(dgNonConformalEdge "dgNonConformalEdge.sh")
set_tests_properties(dgNonConformalEdge
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformalEdge")


#=============================================================================
# dgNonConformalFluids test
#=============================================================================
add_test(dgNonConformalFluids "dgNonConformalFluids.sh")
set_tests_properties(dgNonConformalFluids
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformalFluids")


#=============================================================================
# dgNonConformalFluidsEdge test
#=============================================================================
add_test(dgNonConformalFluidsEdge "dgNonConformalFluidsEdge.sh")
set_tests_properties(dgNonConformalFluidsEdge
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformalFluidsEdge")


#=============================================================================
# dgNonConformal3dFluids test
#=============================================================================
add_test(dgNonConformal3dFluids "dgNonConformal3dFluids.sh")
set_tests_properties(dgNonConformal3dFluids
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformal3dFluids")


#=============================================================================
# dgNonConformal3dFluidsP1P2 test
#=============================================================================
add_test(dgNonConformal3dFluidsP1P2 "dgNonConformal3dFluidsP1P2.sh")
set_tests_properties(dgNonConformal3dFluidsP1P2
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 8
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformal3dFluidsP1P2")


#=============================================================================
# dgNonConformal3dFluidsHexTet test
#=============================================================================
add_test(dgNonConformal3dFluidsHexTet "dgNonConformal3dFluidsHexTet.sh")
set_tests_properties(dgNonConformal3dFluidsHexTet
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformal3dFluidsHexTet")


#=============================================================================
# dgNonConformalThreeBlade test
#=============================================================================
add_test(dgNonConformalThreeBlade "dgNonConformalThreeBlade.sh")
set_tests_properties(dgNonConformalThreeBlade
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformalThreeBlade")


#=============================================================================
# overset test
#=============================================================================
add_test(overset "overset.sh")
set_tests_properties(overset
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 6
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/overset")


#=============================================================================
# oversetFluids test
#=============================================================================
add_test(oversetFluids "oversetFluids.sh")
set_tests_properties(oversetFluids
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 6
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/oversetFluids")


#=============================================================================
# oversetFluidsEdge test
#=============================================================================
add_test(oversetFluidsEdge "oversetFluidsEdge.sh")
set_tests_properties(oversetFluidsEdge
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 6
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/oversetFluidsEdge")


#=============================================================================
# concentricRad test
#=============================================================================
add_test(concentricRad "concentricRad.sh")
set_tests_properties(concentricRad
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/concentricRad")


#=============================================================================
# movingCylinder test
#=============================================================================
add_test(movingCylinder "movingCylinder.sh")
set_tests_properties(movingCylinder
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/movingCylinder")


#=============================================================================
# elemBackStepLRSST test
#=============================================================================
add_test(elemBackStepLRSST "elemBackStepLRSST.sh")
set_tests_properties(elemBackStepLRSST
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/elemBackStepLRSST")


#=============================================================================
# ductWedge test
#=============================================================================
add_test(ductWedge "ductWedge.sh")
set_tests_properties(ductWedge
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 2
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/ductWedge")


#=============================================================================
# heatedBackStep test
#=============================================================================
add_test(heatedBackStep "heatedBackStep.sh")
set_tests_properties(heatedBackStep
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/heatedBackStep")


#=============================================================================
# edgePipeCHT test
#=============================================================================
add_test(edgePipeCHT "edgePipeCHT.sh")
set_tests_properties(edgePipeCHT
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/edgePipeCHT")


#=============================================================================
# elemPipeCHT test
#=============================================================================
add_test(elemPipeCHT "elemPipeCHT.sh")
set_tests_properties(elemPipeCHT
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/elemPipeCHT")


#=============================================================================
# heliumPlume test (with restart; mixed edge/element)
#=============================================================================
add_test(heliumPlume "heliumPlume.sh")
set_tests_properties(heliumPlume
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 8
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/heliumPlume")


#=============================================================================
# edgeContact3D test (cylinder)
#=============================================================================
add_test(edgeContact3D "edgeContact3D.sh")
set_tests_properties(edgeContact3D
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 8
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/edgeContact3D")


#=============================================================================
# dgNonConformalEdgeCylinder test
#=============================================================================
add_test(dgNonConformalEdgeCylinder "dgNonConformalEdgeCylinder.sh")
set_tests_properties(dgNonConformalEdgeCylinder
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 8
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformalEdgeCylinder")


#=============================================================================
# fluidsPmrChtPeriodic test
#=============================================================================
add_test(fluidsPmrChtPeriodic "fluidsPmrChtPeriodic.sh")
set_tests_properties(fluidsPmrChtPeriodic
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 8
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/fluidsPmrChtPeriodic")


#=============================================================================
# nonIsoElemOpenJet test
#=============================================================================
add_test(nonIsoElemOpenJet "nonIsoElemOpenJet.sh")
set_tests_properties(nonIsoElemOpenJet
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/nonIsoElemOpenJet")


#=============================================================================
# nonIsoEdgeOpenJet test
#=============================================================================
add_test(nonIsoEdgeOpenJet "nonIsoEdgeOpenJet.sh")
set_tests_properties(nonIsoEdgeOpenJet
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/nonIsoEdgeOpenJet")


#=============================================================================
# hdf5VarZChi test
#=============================================================================
# add_test(hdf5VarZChi "hdf5VarZChi.sh")
# set_tests_properties(hdf5VarZChi
#   PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
#               FAIL_REGULAR_EXPRESSION "FAILED"
#               TIMEOUT 400
#               PROCESSORS 4
#   WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/hdf5VarZChi")


#=============================================================================
# elemHybridFluids test
#=============================================================================
add_test(elemHybridFluids "elemHybridFluids.sh")
set_tests_properties(elemHybridFluids
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 8
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/elemHybridFluids")


#=============================================================================
# elemHybridFluidsShift test
#=============================================================================
add_test(elemHybridFluidsShift "elemHybridFluidsShift.sh")
set_tests_properties(elemHybridFluidsShift
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 8
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/elemHybridFluidsShift")


#=============================================================================
# edgeHybridFluids test
#=============================================================================
add_test(edgeHybridFluids "edgeHybridFluids.sh")
set_tests_properties(edgeHybridFluids
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 8
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/edgeHybridFluids")


#=============================================================================
# elemClosedDomain test
#=============================================================================
add_test(elemClosedDomain "elemClosedDomain.sh")
set_tests_properties(elemClosedDomain
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 2
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/elemClosedDomain")


#=============================================================================
# mixedTetPipe test
#=============================================================================
add_test(mixedTetPipe "mixedTetPipe.sh")
set_tests_properties(mixedTetPipe
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 8
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/mixedTetPipe")


#=============================================================================
# inputFireEdgeUpwind test
#=============================================================================
add_test(inputFireEdgeUpwind "inputFireEdgeUpwind.sh")
set_tests_properties(inputFireEdgeUpwind
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/inputFireEdgeUpwind")


#=============================================================================
# inputFireElem test
#=============================================================================
add_test(inputFireElem "inputFireElem.sh")
set_tests_properties(inputFireElem
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/inputFireElem")


#=============================================================================
# nonIsoNonUniformElemOpenJet test
#=============================================================================
add_test(nonIsoNonUniformElemOpenJet "nonIsoNonUniformElemOpenJet.sh")
set_tests_properties(nonIsoNonUniformElemOpenJet
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/nonIsoNonUniformElemOpenJet")


#=============================================================================
# nonIsoNonUniformEdgeOpenJet test
#=============================================================================
add_test(nonIsoNonUniformEdgeOpenJet "nonIsoNonUniformEdgeOpenJet.sh")
set_tests_properties(nonIsoNonUniformEdgeOpenJet
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/nonIsoNonUniformEdgeOpenJet")


#=============================================================================
# milestoneRun test
#=============================================================================
add_test(milestoneRun "milestoneRun.sh")
set_tests_properties(milestoneRun
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/milestoneRun")


#=============================================================================
# heatedWaterChannel test
#=============================================================================
add_test(heatedWaterChannel "heatedWaterChannel.sh")
set_tests_properties(heatedWaterChannel
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/heatedWaterChannel")


#=============================================================================
# variableDensMMS test (edge and element)
#=============================================================================
add_test(variableDensMMS "variableDensMMS.sh")
set_tests_properties(variableDensMMS
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 2
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/variableDensMMS")


#=============================================================================
# actuatorLine test
#=============================================================================
add_test(actuatorLine "actuatorLine.sh")
set_tests_properties(actuatorLine
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 8
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/actuatorLine")


#=============================================================================
# femHC test
#=============================================================================
add_test(femHC "femHC.sh")
set_tests_properties(femHC
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 2
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/femHC")


#=============================================================================
# ablUnstableEdge test
#=============================================================================
add_test(ablUnstableEdge "ablUnstableEdge.sh")
set_tests_properties(ablUnstableEdge
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/ablUnstableEdge")


#=============================================================================
# ablStableElem test
#=============================================================================
add_test(ablStableElem "ablStableElem.sh")
set_tests_properties(ablStableElem
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 4
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/ablStableElem")


#=============================================================================
# unit tests
#=============================================================================
add_test(unitTests "run_unit_tests.sh")
set_tests_properties(unitTests
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 400
              PROCESSORS 2
  WORKING_DIRECTORY "${NIGHTLY_TEST_RESULT_DIRECTORY}/unitTests")

#=============================================================================
#
# Performance tests
#
#=============================================================================

#=============================================================================
# waleElemXflowMixFrac3.5m test
#=============================================================================
add_test(waleElemXflowMixFrac3.5m "waleElemXflowMixFrac3.5m.sh")
set_tests_properties(waleElemXflowMixFrac3.5m
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 600
              PROCESSORS 8
  WORKING_DIRECTORY "${PERF_TEST_RESULT_DIRECTORY}/waleElemXflowMixFrac3.5m")


#=============================================================================
# uqSlidingMesh test
#=============================================================================
add_test(uqSlidingMesh "uqSlidingMesh.sh")
set_tests_properties(uqSlidingMesh
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 1000
              PROCESSORS 8
  WORKING_DIRECTORY "${PERF_TEST_RESULT_DIRECTORY}/uqSlidingMesh")


#=============================================================================
# uqSlidingMeshDG test
#=============================================================================
add_test(uqSlidingMeshDG "uqSlidingMeshDG.sh")
set_tests_properties(uqSlidingMeshDG
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 1000
              PROCESSORS 8
  WORKING_DIRECTORY "${PERF_TEST_RESULT_DIRECTORY}/uqSlidingMeshDG")


#=============================================================================
# oversetHybrid test
#=============================================================================
add_test(oversetHybrid "oversetHybrid.sh")
set_tests_properties(oversetHybrid
  PROPERTIES  PASS_REGULAR_EXPRESSION "PASSED"
              FAIL_REGULAR_EXPRESSION "FAILED"
              TIMEOUT 1000
              PROCESSORS 8
  WORKING_DIRECTORY "${PERF_TEST_RESULT_DIRECTORY}/oversetHybrid")
