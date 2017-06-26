/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "UnitTestRealm.h"
#include "UnitTestUtils.h"

#include "LinearSolvers.h"
#include "Realms.h"
#include "Realm.h"
#include "InputOutputRealm.h"
#include "SolutionOptions.h"
#include "TimeIntegrator.h"
#include "ComputeSSTMaxLengthScaleElemAlgorithm.h"

#include <string>

namespace {

const std::string naluDefaultInputs =
  "Simulations:                                                            \n"
  "  - name: sim1                                                          \n"
  "    time_integrator: ti_1                                               \n"
  "    optimizer: opt1                                                     \n"
  "                                                                        \n"
  "linear_solvers:                                                         \n"
  "                                                                        \n"
  "  - name: solve_scalar                                                  \n"
  "    type: tpetra                                                        \n"
  "    method: gmres                                                       \n"
  "    preconditioner: sgs                                                 \n"
  "    tolerance: 1e-5                                                     \n"
  "    max_iterations: 50                                                  \n"
  "    kspace: 50                                                          \n"
  "    output_level: 0                                                     \n"
  "                                                                        \n"
  "  - name: solve_cont                                                    \n"
  "    type: tpetra                                                        \n"
  "    method: gmres                                                       \n"
  "    preconditioner: muelu                                               \n"
  "    tolerance: 1e-5                                                     \n"
  "    max_iterations: 50                                                  \n"
  "    kspace: 50                                                          \n"
  "    output_level: 0                                                     \n"
  "    recompute_preconditioner: no                                        \n"
  "    muelu_xml_file_name: milestone.xml                                  \n"
  "                                                                        \n"
  "Time_Integrators:                                                       \n"
  "  - StandardTimeIntegrator:                                             \n"
  "      name: ti_1                                                        \n"
  "      start_time: 0                                                     \n"
  "      time_step: 0.1                                                    \n"
  "      termination_time: 1.5                                             \n"
  "      time_stepping_type: adaptive                                      \n"
  "      time_step_count: 0                                                \n"
  "      second_order_accuracy: yes                                        \n"
  "                                                                        \n"
  "      realms: []                                                        \n"
  "                                                                        \n"
  ;

const std::string realmDefaultSettings =
  "- name: unitTestRealm                                                  \n"
  "  use_edges: no                                                        \n"
  "                                                                       \n"
  "  time_step_control:                                                   \n"
  "    target_courant: 2.0                                                \n"
  "    time_step_change_factor: 1.2                                       \n"
  "                                                                       \n"
  "  solution_options:                                                    \n"
  "    name: unitTestRealmOptions                                         \n"
  "    turbulence_model: laminar                                          \n"
  "    interp_rhou_together_for_mdot: yes                                 \n"
  "    use_consolidated_solver_algorithm: no                              \n"
  "    reduced_sens_cvfem_poisson: no                                     \n"
  "                                                                       \n"
  "    options:                                                           \n"
  "      - laminar_prandtl:                                               \n"
  "          enthalpy: 0.7                                                \n"
  "      - turbulent_prandtl:                                             \n"
  "          enthalpy: 1.0                                                \n"
  "      - shifted_gradient_operator:                                     \n"
  "          velocity: no                                                 \n"
  "          pressure: no                                                 \n"
  "          mixture_fraction: no                                         \n"
  ;
}

namespace unit_test_utils {

YAML::Node get_default_inputs() {
  YAML::Node doc = YAML::Load(naluDefaultInputs);

  return doc;
}

YAML::Node get_realm_default_node() {
  YAML::Node node = YAML::Load(realmDefaultSettings);
  return node[0];
}

NaluTest::NaluTest(YAML::Node& doc)
  : doc_(doc),
    comm_(MPI_COMM_WORLD),
    spatialDim_(3),
    sim_(doc_)
{
  // NaluEnv log file
  std::string logFileName = "unittestX_naluwrapper.log";
  auto testInfo = ::testing::UnitTest::GetInstance()->current_test_info();
  if (testInfo) {
    std::string caseName = testInfo->test_case_name();
    std::string caseInstance = testInfo->name();
    logFileName = caseName + "." + caseInstance + ".log";
  }
  sierra::nalu::NaluEnv::self().set_log_file_stream(logFileName, false);

  sim_.linearSolvers_ = new sierra::nalu::LinearSolvers(sim_);
  sim_.realms_ = new sierra::nalu::Realms(sim_);
  sim_.timeIntegrator_ = new sierra::nalu::TimeIntegrator(&sim_);

  sim_.linearSolvers_->load(doc_);
  sim_.timeIntegrator_->load(doc_);
}

sierra::nalu::Realm&
NaluTest::create_realm(const YAML::Node& realm_node, const std::string realm_type)
{
  sierra::nalu::Realm* realm = nullptr;
  if (realm_type == "multi_physics")
    realm = new sierra::nalu::Realm(*sim_.realms_, realm_node);
  else
    realm = new sierra::nalu::InputOutputRealm(*sim_.realms_, realm_node);

  // Populate solution options
  realm->solutionOptions_->load(realm_node);

  // Set-up mesh metadata and bulkdata ... let user fill mesh in test function
  realm->metaData_ = new stk::mesh::MetaData(spatialDim_);
  realm->bulkData_ = new stk::mesh::BulkData(*realm->metaData_, comm_);
  sim_.realms_->realmVector_.push_back(realm);

  return *realm;
}

TEST(NaluMock, test_nalu_mock)
{
  // 1. Create dummy simulation object
  YAML::Node doc = get_default_inputs();
  unit_test_utils::NaluTest naluObj(doc);

  // 2. Create a Realm input node used to fill data
  const YAML::Node realm_node = get_realm_default_node();
  // Modify the default node for the test

  // 3. Create the Realm node
  sierra::nalu::Realm& realm = naluObj.create_realm(realm_node);

  // 4. Create mesh and get the default part for registration with Algorithm
  unit_test_utils::fill_hex8_mesh("generated:10x10x10", *realm.bulkData_);
  stk::mesh::Part* meshPart = realm.metaData_->get_part("block_1");

  // 5. Create necessary fields...

  // 6. Initialize the Algorithm to be tested...
  sierra::nalu::ComputeSSTMaxLengthScaleElemAlgorithm sstAlg(realm, meshPart);

  // 7. Perform tests

}

}  // unit_test_utils
