#include <gtest/gtest.h>

#include <stk_util/environment/WallTime.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include "UnitTestUtils.h"
#include "UnitTestHelperObjects.h"

#include <SolutionOptions.h>
#include <TimeIntegrator.h>
#include <SolverAlgorithm.h>
#include <SolverAlgorithmDriver.h>
#include <MomentumAdvDiffElemKernel.h>
#include <AssembleContinuityElemOpenSolverAlgorithm.h>
#include <AssembleMomentumElemOpenSolverAlgorithm.h>
#include <ElemDataRequests.h>

class AssembleAlgorithmMesh : public Hex8MeshWithNSOFields {};

void declare_open_solver_alg_fields(stk::mesh::MetaData& meta)
{
  GenericFieldType& openMassFlowRateField = meta.declare_field<GenericFieldType>(meta.side_rank(), "open_mass_flow_rate");
  GenericFieldType& exposedAreaVecField = meta.declare_field<GenericFieldType>(meta.side_rank(), "exposed_area_vector");
  GenericFieldType& dudxField = meta.declare_field<GenericFieldType>(stk::topology::NODE_RANK, "dudx");
  VectorFieldType& velBcField = meta.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "open_velocity_bc");
  ScalarFieldType& pressBcField = meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure_bc");
  double zero = 0.0;
  stk::mesh::put_field(exposedAreaVecField, meta.universal_part(), 1, &zero);
  stk::mesh::put_field(openMassFlowRateField, meta.universal_part(), 1, &zero);
  stk::mesh::put_field(velBcField, meta.universal_part(), 1, &zero);
  stk::mesh::put_field(dudxField, meta.universal_part(), 1, &zero);
  stk::mesh::put_field(pressBcField, meta.universal_part(), 1, &zero);
}

TEST_F(AssembleAlgorithmMesh, onlyAssembleRhs)
{
  declare_open_solver_alg_fields(meta);
  fill_mesh_and_initialize_test_fields("generated:2x2x2|sideset:xXyYzZ");

  sierra::nalu::SolutionOptions solnOpts;
  solnOpts.meshMotion_ = false;
  solnOpts.meshDeformation_ = false;
  solnOpts.externalMeshDeformation_ = false;
  solnOpts.includeDivU_ = 0.0;

  const unsigned numDof = 3;
  stk::mesh::Part* block_1 = meta.get_part("block_1");
  unit_test_utils::HelperObjectsNewME helperObjs(bulk, stk::topology::HEX_8, numDof, block_1);

  sierra::nalu::TimeIntegrator timeIntegrator;
  timeIntegrator.timeStepN_ = 0.1;
  timeIntegrator.timeStepNm1_ = 0.1;
  timeIntegrator.gamma1_ = 1.0;
  timeIntegrator.gamma2_ = -1.0;
  timeIntegrator.gamma3_ = 0.0;

  helperObjs.realm.timeIntegrator_ = &timeIntegrator;

  std::unique_ptr<sierra::nalu::Kernel> advDiffKernel(
     new sierra::nalu::MomentumAdvDiffElemKernel<sierra::nalu::AlgTraitsHex8>(
        bulk, solnOpts, velocity, viscosity,
        helperObjs.assembleElemSolverAlg->dataNeededByKernels_));

  helperObjs.assembleElemSolverAlg->activeKernels_.push_back(advDiffKernel.get());

  sierra::nalu::AssembleMomentumElemOpenSolverAlgorithm assembleMomentumElemOpenSolverAlgorithm(helperObjs.realm, block_1, &helperObjs.eqSystem);

  sierra::nalu::AssembleContinuityElemOpenSolverAlgorithm assembleContinuityElemOpenSolverAlgorithm(helperObjs.realm, block_1, &helperObjs.eqSystem);

  bool onlyRhs = true;

  helperObjs.assembleElemSolverAlg->set_only_assemble_rhs(onlyRhs);
  helperObjs.assembleElemSolverAlg->execute();

  EXPECT_TRUE(helperObjs.linsys->ignoredLhs_);
  helperObjs.linsys->ignoredLhs_ = false;

  assembleMomentumElemOpenSolverAlgorithm.set_only_assemble_rhs(onlyRhs);
  assembleMomentumElemOpenSolverAlgorithm.execute();

  EXPECT_TRUE(helperObjs.linsys->ignoredLhs_);
  helperObjs.linsys->ignoredLhs_ = false;

  assembleContinuityElemOpenSolverAlgorithm.set_only_assemble_rhs(onlyRhs);
  assembleContinuityElemOpenSolverAlgorithm.execute();

  EXPECT_TRUE(helperObjs.linsys->ignoredLhs_);
}

class TestSolverAlg : public sierra::nalu::SolverAlgorithm {
  public:
    TestSolverAlg(sierra::nalu::Realm& realm, stk::mesh::Part* part,
                  sierra::nalu::EquationSystem* eqSystem)
     : SolverAlgorithm(realm, part, eqSystem)
    {}

    void execute() {}
    void initialize_connectivity() {}
};

TEST_F(AssembleAlgorithmMesh, SolverAlgDriver_onlyAssembleRhs)
{
  declare_open_solver_alg_fields(meta);
  fill_mesh_and_initialize_test_fields("generated:2x2x2|sideset:xXyYzZ");

  const unsigned numDof = 3;
  stk::mesh::Part* block_1 = meta.get_part("block_1");
  unit_test_utils::HelperObjectsNewME helperObjs(bulk, stk::topology::HEX_8, numDof, block_1);

  //these algs will be deleted by SolverAlgorithmDriver
  TestSolverAlg* alg1 = new TestSolverAlg(helperObjs.realm, block_1, &helperObjs.eqSystem);
  TestSolverAlg* alg2 = new TestSolverAlg(helperObjs.realm, block_1, &helperObjs.eqSystem);
  TestSolverAlg* alg3 = new TestSolverAlg(helperObjs.realm, block_1, &helperObjs.eqSystem);
  TestSolverAlg* alg4 = new TestSolverAlg(helperObjs.realm, block_1, &helperObjs.eqSystem);

  bool onlyRhs = true;

  sierra::nalu::SolverAlgorithmDriver solverAlgDriver(helperObjs.realm);

  solverAlgDriver.solverAlgorithmMap_.insert(std::make_pair("alg1", alg1));
  solverAlgDriver.solverAlgMap_.insert(std::make_pair(sierra::nalu::INFLOW, alg2));
  solverAlgDriver.solverConstraintAlgMap_.insert(std::make_pair(sierra::nalu::INFLOW, alg3));
  solverAlgDriver.solverDirichAlgMap_.insert(std::make_pair(sierra::nalu::INFLOW, alg4));

  solverAlgDriver.set_only_assemble_rhs(onlyRhs);

  EXPECT_EQ(onlyRhs, alg1->get_only_assemble_rhs());
  EXPECT_EQ(onlyRhs, alg2->get_only_assemble_rhs());
  EXPECT_EQ(onlyRhs, alg3->get_only_assemble_rhs());
  EXPECT_EQ(onlyRhs, alg4->get_only_assemble_rhs());
}

