# Prepare the tests directories

#=============================================================================
#
# Prepare test directory
#
#=============================================================================

# make a clean test directory
file(REMOVE_RECURSE ${TEST_RESULT_DIRECTORY})
file(MAKE_DIRECTORY ${TEST_RESULT_DIRECTORY})

# copy executables to the test directory
file(COPY ${CTEST_NALU_BINARY_NAME}
          ${CTEST_UNITTEST_BINARY_NAME}
          ${NALURTEST_DIR}/pass_fail.sh
     DESTINATION ${TEST_RESULT_DIRECTORY})


#=============================================================================
#
# Nightly directories
#
#=============================================================================

#=============================================================================
# periodic3dElem test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/periodic3dElem"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_MESH_DIR}/periodic3d.g" 
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/periodic3dElem")


#=============================================================================
# periodic3dEdge test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/periodic3dEdge"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_MESH_DIR}/periodic3d.g" 
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/periodic3dEdge")


#=============================================================================
# quad9HC test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/quad9HC"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/100x50_P2n.g.*")
file(COPY               ${meshes}
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/quad9HC")


#=============================================================================
# steadyTaylorVortex test
#=============================================================================i
file(COPY               "${NALURTEST_NIGHTLY_DIR}/steadyTaylorVortex"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/steadyTaylorVortex")


#=============================================================================
# hoVortex test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/hoVortex"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/100x50_P2n.g.*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/hoVortex")

   
#=============================================================================
# hoHelium test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/hoHelium"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/hoHelium")


#=============================================================================
# dgNonConformal test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/dgNonConformal"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/2DTwoBlock_R0_rev.g.*")
file(COPY               ${meshes}
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformal")


#=============================================================================
# dgNonConformalEdge test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/dgNonConformalEdge"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/2DTwoBlock_R0_rev.g.*")
file(COPY               ${meshes}
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformalEdge")


#=============================================================================
# dgNonConformalFluids test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/dgNonConformalFluids"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/NACA.g.*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformalFluids")


#=============================================================================
# dgNonConformalFluidsEdge test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/dgNonConformalFluidsEdge"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/NACA.g.*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformalFluidsEdge")


#=============================================================================
# dgNonConformal3dFluids test
#=============================================================================i
file(COPY               "${NALURTEST_NIGHTLY_DIR}/dgNonConformal3dFluids"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/twoBlockMesh_cgs.g.*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformal3dFluids")


#=============================================================================
# dgNonConformal3dFluidsP1P2 test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/dgNonConformal3dFluidsP1P2"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformal3dFluidsP1P2")


#=============================================================================
# dgNonConformal3dFluidsHexTet test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/dgNonConformal3dFluidsHexTet"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/twoBlockMeshHexTet_cgs.g.*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformal3dFluidsHexTet")


#=============================================================================
# dgNonConformalThreeBlade test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/dgNonConformalThreeBlade"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_MESH_DIR}/threeBladeMesh.g"
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformalThreeBlade")


#=============================================================================
# overset test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/overset"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/oversetMeshAligned.g*")
file(COPY               ${meshes}
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/overset")


#=============================================================================
# oversetFluids test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/oversetFluids"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/oversetMeshAligned.g*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/oversetFluids")


#=============================================================================
# oversetFluidsEdge test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/oversetFluidsEdge"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/oversetMeshAligned.g*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/oversetFluidsEdge")


#=============================================================================
# concentricRad test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/concentricRad"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")


#=============================================================================
# movingCylinder test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/movingCylinder"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_XML_DIR}/milestone_aspect_ratio.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/movingCylinder")


#=============================================================================
# elemBackStepLRSST test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/elemBackStepLRSST"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/elemBackStepLRSST")


#=============================================================================
# ductWedge test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/ductWedge"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/ductwedge.g.*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/ductWedge")


#=============================================================================
# heatedBackStep test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/heatedBackStep"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/heatedBackStep")


#=============================================================================
# edgePipeCHT test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/edgePipeCHT"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes1       "${NALURTEST_MESH_DIR}/elbow.g.*")
file(GLOB meshes2       "${NALURTEST_MESH_DIR}/horseshoe.g.*")
file(COPY               ${meshes1}
                        ${meshes2}
                        "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/edgePipeCHT")


#=============================================================================
# elemPipeCHT test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/elemPipeCHT"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes1       "${NALURTEST_MESH_DIR}/elbow.g.*")
file(GLOB meshes2       "${NALURTEST_MESH_DIR}/horseshoe.g.*")
file(COPY               ${meshes1}
                        ${meshes2}
                        "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/elemPipeCHT")


#=============================================================================
# heliumPlume test (with restart; mixed edge/element)
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/heliumPlume"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/heliumPlume")


#=============================================================================
# dgNonConformalEdgeCylinder test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/dgNonConformalEdgeCylinder"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/rot_cyl_14.exo*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/dgNonConformalEdgeCylinder")


#=============================================================================
# fluidsPmrChtPeriodic test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/fluidsPmrChtPeriodic"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/fluidsPmrChtPeriodic")


#=============================================================================
# nonIsoElemOpenJet test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/nonIsoElemOpenJet"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/2cm_ped_35K_mks.g*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/nonIsoElemOpenJet")


#=============================================================================
# nonIsoEdgeOpenJet test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/nonIsoEdgeOpenJet"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/2cm_ped_35K_mks.g*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/nonIsoEdgeOpenJet")


#=============================================================================
# hdf5VarZChi test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/hdf5VarZChi"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/2cm_ped_35K_mks.g*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/hdf5VarZChi")


#=============================================================================
# elemHybridFluids test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/elemHybridFluids"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/hybrid.g*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/milestone_aspect_ratio.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/elemHybridFluids")


#=============================================================================
# elemHybridFluidsShift test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/elemHybridFluidsShift"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/hybrid.g*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/elemHybridFluidsShift")


#=============================================================================
# edgeHybridFluids test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/edgeHybridFluids"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/hybrid.g*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/edgeHybridFluids")


#=============================================================================
# elemClosedDomain test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/elemClosedDomain"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/elemClosedDomain")


#=============================================================================
# mixedTetPipe test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/mixedTetPipe"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/pipeTet.g*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/mixedTetPipe")


#=============================================================================
# inputFireEdgeUpwind test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/inputFireEdgeUpwind"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/13K_oneWay.e.*")
file(COPY               ${meshes}
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/inputFireEdgeUpwind")


#=============================================================================
# inputFireElem test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/inputFireElem"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/13K_oneWay.e.*")
file(COPY               ${meshes}
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/inputFireElem")


#=============================================================================
# nonIsoNonUniformElemOpenJet test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/nonIsoNonUniformElemOpenJet"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/2cm_ped_35K_mks.g*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/nonIsoNonUniformElemOpenJet")


#=============================================================================
# nonIsoNonUniformEdgeOpenJet test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/nonIsoNonUniformEdgeOpenJet"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/2cm_ped_35K_mks.g*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/matches_ml_default.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/nonIsoNonUniformEdgeOpenJet")


#=============================================================================
# milestoneRun test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/milestoneRun"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_MESH_DIR}/1cm_ped_35K.g"
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/milestoneRun")


#=============================================================================
# heatedWaterChannel test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/heatedWaterChannel"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_XML_DIR}/milestone.xml"
                        "${NALURTEST_XML_DIR}/milestone_aspect_ratio.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/heatedWaterChannel")


#=============================================================================
# variableDensMMS test (edge and element)
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/variableDensMMS"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/variableDensMMS")


#=============================================================================
# actuatorLine test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/actuatorLine"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/actuatorLine")


#=============================================================================
# femHC test
#=============================================================================i
file(COPY               "${NALURTEST_NIGHTLY_DIR}/femHC"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_MESH_DIR}/periodic3d.g"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/femHC")


#=============================================================================
# cvfemHC test
#=============================================================================i
file(COPY               "${NALURTEST_NIGHTLY_DIR}/cvfemHC"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_MESH_DIR}/rot_cyl_14.exo*"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/cvfemHC")


#=============================================================================
# ablUnstableEdge test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/ablUnstableEdge"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_MESH_DIR}/abl_1km_cube_toy.g"
                        "${NALURTEST_MESH_DIR}/abl_io.g"
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/ablUnstableEdge")


#=============================================================================
# ablStableElem test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/ablStableElem"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_MESH_DIR}/abl_1km_cube_toy.g"
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/ablStableElem")


#=============================================================================
# ekmanSpiral test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/ekmanSpiral"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_MESH_DIR}/ekmanSpiral.g"
                        "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}/ekmanSpiral")

#=============================================================================
# dgMMS test
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/dgMMS"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")

#=============================================================================
# unit tests
#=============================================================================
file(COPY               "${NALURTEST_NIGHTLY_DIR}/unitTests"
     DESTINATION        "${NIGHTLY_TEST_RESULT_DIRECTORY}")


#=============================================================================
#
# Performance tests
#
#=============================================================================

#=============================================================================
# waleElemXflowMixFrac3.5m test
#=============================================================================
file(COPY               "${NALURTEST_PERF_DIR}/waleElemXflowMixFrac3.5m"
     DESTINATION        "${PERF_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_XML_DIR}/milestone_aspect_ratio_smooth.xml"
     DESTINATION        "${PERF_TEST_RESULT_DIRECTORY}/waleElemXflowMixFrac3.5m")


#=============================================================================
# uqSlidingMeshDG test
#=============================================================================
file(COPY               "${NALURTEST_PERF_DIR}/uqSlidingMeshDG"
     DESTINATION        "${PERF_TEST_RESULT_DIRECTORY}")

file(GLOB meshes        "${NALURTEST_MESH_DIR}/uqvawt_corrected.exo.*")
file(COPY               ${meshes}
                        "${NALURTEST_XML_DIR}/milestone_aspect_ratio_smooth.xml"
     DESTINATION        "${PERF_TEST_RESULT_DIRECTORY}/uqSlidingMeshDG")


#=============================================================================
# oversetHybrid test
#=============================================================================
file(COPY               "${NALURTEST_PERF_DIR}/oversetHybrid"
     DESTINATION        "${PERF_TEST_RESULT_DIRECTORY}")

file(COPY               "${NALURTEST_XML_DIR}/milestone.xml"
     DESTINATION        "${PERF_TEST_RESULT_DIRECTORY}/oversetHybrid")


