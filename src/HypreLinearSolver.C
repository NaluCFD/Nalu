/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifdef NALU_USES_HYPRE

#include "HypreLinearSolver.h"
#include <stk_util/environment/ReportHandler.hpp>

#include "NaluHypre.h"
#include "PeriodicManager.h"
#include "NonConformalManager.h"
#include "overset/OversetManager.h"

#include "FieldTypeDef.h"

#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/Field.hpp"

namespace sierra {
namespace nalu {

HypreLinearSolver::HypreLinearSolver(
  std::string name,
  HypreLinearSolverConfig* config,
  LinearSolvers* linearSolvers)
  : LinearSolver(name, linearSolvers, config)
{}

HypreLinearSolver::~HypreLinearSolver()
{
  destroyLinearSolver();
}

void
HypreLinearSolver::destroyLinearSolver()
{
  // if (!hypreSolver_.is_null())
  //   hypreSolver_->Destroy();
  hypreSolver_ = Teuchos::null;
}

void
HypreLinearSolver::setupLinearSolver(
  Teuchos::RCP<LinSys::Vector> /* sln */,
  Teuchos::RCP<LinSys::Matrix> matrix,
  Teuchos::RCP<LinSys::Vector> rhs,
  Teuchos::RCP<LinSys::MultiVector>)
{
  ThrowRequire(!matrix.is_null());
  ThrowRequire(!rhs.is_null());

  matrix_ = matrix;
  rhs_ = rhs;

}

int
HypreLinearSolver::solve(
  Teuchos::RCP<LinSys::Vector> sln,
  int& iters,
  double& finalResidNrm,
  bool)
{
  ThrowRequire(!sln.is_null());

  const int status = 0;
  iters = 0;
  finalResidNrm = 0.0;

  if (hypreSolver_.is_null()) {
    sync_hypre_id();
    // Create the solver
    hypreSolver_ = Teuchos::rcp(
      new Ifpack2::NaluHypre<LinSys::Scalar, LinSys::LocalOrdinal, LinSys::GlobalOrdinal, LinSys::Node>(matrix_, realm_, numDof_));

    hypreSolver_->dumpMatrixFiles(config_->getWriteMatrixFiles());
    hypreSolver_->setParameters(*(config_->paramsPrecond()));
    hypreSolver_->initialize();
    hypreSolver_->compute();

    if (realm_->bulk_data().parallel_rank() == 0)
      std::cerr << "Intializing HYPRE Solver" << std::endl;
  }

  hypreSolver_->apply(*rhs_, *sln);

  timerPrecond_ = (hypreSolver_->getInitializeTime() + hypreSolver_->getComputeTime());

  // auto slvtmp = reinterpret_cast<Ifpack2::NaluHypre<LinSys::Scalar, LinSys::LocalOrdinal, LinSys::GlobalOrdinal, LinSys::Node>*>(hypreSolver_.get());
  iters = hypreSolver_->hypreNumLinIters_;
  finalResidNrm = hypreSolver_->hypreFinalResidual_;

  return status;
}

void
HypreLinearSolver::sync_hypre_id()
{
  auto& bulk = realm_->bulk_data();
  std::vector<const stk::mesh::FieldBase*> fVec{realm_->hypreGlobalId_};

  stk::mesh::copy_owned_to_shared(bulk, fVec);
  stk::mesh::communicate_field_data(bulk.aura_ghosting(), fVec);

  if (realm_->oversetManager_ != nullptr &&
      realm_->oversetManager_->oversetGhosting_ != nullptr)
    stk::mesh::communicate_field_data(
      *realm_->oversetManager_->oversetGhosting_, fVec);

  if (realm_->nonConformalManager_ != nullptr &&
      realm_->nonConformalManager_->nonConformalGhosting_ != nullptr)
    stk::mesh::communicate_field_data(
      *realm_->nonConformalManager_->nonConformalGhosting_, fVec);

  if (realm_->periodicManager_ != nullptr &&
      realm_->periodicManager_->periodicGhosting_ != nullptr)
    realm_->periodicManager_->periodic_parallel_communicate_field(
      realm_->hypreGlobalId_);
}

}  // nalu
}  // sierra

#endif
