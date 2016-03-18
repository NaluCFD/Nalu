/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <LinearSolver.h>
#include <LinearSolvers.h>

#include <NaluEnv.h>
#include <LinearSolverTypes.h>

#include <stk_util/environment/ReportHandler.hpp>

#include <Epetra_FECrsMatrix.h>
#include <Epetra_FEVector.h>
#include <Epetra_Vector.h>
#include <AztecOO.h>
#include <ml_MultiLevelPreconditioner.h>
#include <ml_include.h>

// Tpetra support
#include <BelosLinearProblem.hpp>
#include <BelosMultiVecTraits.hpp>
#include <BelosOperatorTraits.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosSolverManager.hpp>
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>

#include <Ifpack2_Factory.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_Serial.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>

#include <Teuchos_ParameterXMLFileReader.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_CreateEpetraPreconditioner.hpp>

// stk_util
#include <stk_util/environment/CPUTime.hpp>

#include <iostream>

namespace sierra{
namespace nalu{

Simulation *LinearSolver::root() { return linearSolvers_->root(); }
LinearSolvers *LinearSolver::parent() { return linearSolvers_; }

EpetraLinearSolver::EpetraLinearSolver(
  std::string solverName,
  EpetraLinearSolverConfig *config,
  const int * options,
  const double * params, bool mlFlag,
  bool mueLuFlag,
  const Teuchos::RCP<Teuchos::ParameterList> mlParams,
  LinearSolvers *linearSolvers)
  : LinearSolver(solverName, linearSolvers, config->recomputePreconditioner(), config->reusePreconditioner()),
    config_(config),
    solver_(new AztecOO),
    activateML_(mlFlag),
    activateMueLu_(mueLuFlag),
    mlPreconditioner_(0),
    mlParams_(mlParams)
{
  solver_->SetAllAztecOptions(options);
  solver_->SetAllAztecParams(params);
  solver_->SetOutputStream(NaluEnv::self().naluOutputP0());
  solver_->SetErrorStream(NaluEnv::self().naluOutputP0());
}

EpetraLinearSolver::~EpetraLinearSolver () {
  delete solver_;
}

void
EpetraLinearSolver::setSystemObjects(
  Epetra_FECrsMatrix * matrix,
  Epetra_FEVector * rhs)
{
  ThrowRequire(matrix);
  ThrowRequire(rhs);
  ThrowRequire(solver_);
  solver_->SetUserMatrix(matrix);
  solver_->SetRHS(rhs);
}

void
EpetraLinearSolver::populateMLCoordinates(double * xcoords, double * ycoords, double * zcoords)
{
  if(!activateML_) return;

  mlParams_->set("x-coordinates", xcoords);
  mlParams_->set("y-coordinates", ycoords);
  mlParams_->set("z-coordinates", zcoords);

}
void
EpetraLinearSolver::setMueLuCoordinates(Teuchos::RCP<Epetra_MultiVector> coords)
{
  if(!activateMueLu_) return;

  mueLuCoordinates_ = coords;
}

void
EpetraLinearSolver::destroyML()
{
  if (!activateML_) return;
  delete mlPreconditioner_;
}

void
EpetraLinearSolver::destroyMueLu()
{
  if (!activateMueLu_) return;
  mueLuPreconditioner_ = Teuchos::null;
}

int
EpetraLinearSolver::solve(
  Epetra_Vector * sln,
  int & iteration_count,
  double & scaledResidual)
{
  ThrowRequire(solver_);
  ThrowRequire(sln);
  Epetra_RowMatrix * matrix = solver_->GetUserMatrix();
  ThrowRequire(matrix);
  ThrowRequire(solver_->GetRHS());
  solver_->SetLHS(sln);

  double time = -stk::cpu_time();
  if (activateML_)
  {
    if (mlPreconditioner_ == 0)
      mlPreconditioner_ = new ML_Epetra::MultiLevelPreconditioner(*matrix, *mlParams_);
    solver_->SetPrecOperator(mlPreconditioner_);
  }
  if (activateMueLu_)
  {
    Teuchos::RCP<Teuchos::Time> tm = Teuchos::TimeMonitor::getNewTimer("nalu MueLu preconditioner setup");
    Teuchos::TimeMonitor timeMon(*tm);

    if (recomputePreconditioner_ || mueLuPreconditioner_ == Teuchos::null)
    {
      std::string xmlFileName = config_->muelu_xml_file();
      mueLuPreconditioner_ = MueLu::CreateEpetraPreconditioner(mueLuMat_, xmlFileName, mueLuCoordinates_);
    }
    else if (reusePreconditioner_) {
      MueLu::ReuseEpetraPreconditioner(mueLuMat_, *mueLuPreconditioner_);
    }
    if (config_->getSummarizeMueluTimer())
      Teuchos::TimeMonitor::summarize(std::cout, false, true, false, Teuchos::Union);

    solver_->SetPrecOperator(mueLuPreconditioner_.getRawPtr());
  }
  time += stk::cpu_time();
  timerPrecond_ += time;
  const int max_iterations = solver_->GetAztecOption(AZ_max_iter);
  const double tol = solver_->GetAllAztecParams()[AZ_tol];

  const int status = solver_->Iterate(max_iterations, tol);
  iteration_count = solver_->NumIters();
  scaledResidual = solver_->ScaledResidual();

  if (activateML_ && recomputePreconditioner_)
  {
    delete mlPreconditioner_; mlPreconditioner_ = 0;
  }

  return status;
}

TpetraLinearSolver::TpetraLinearSolver(
  std::string solverName,
  TpetraLinearSolverConfig *config,
  const Teuchos::RCP<Teuchos::ParameterList> params,
  const Teuchos::RCP<Teuchos::ParameterList> paramsPrecond,
  LinearSolvers *linearSolvers)
  : LinearSolver(solverName,linearSolvers, config->recomputePreconditioner(), config->reusePreconditioner()),
    config_(config),
    params_(params),
    paramsPrecond_(paramsPrecond),
    activateMueLu_(config->use_MueLu())
{
}

TpetraLinearSolver::~TpetraLinearSolver()
{
  destroyLinearSolver();
}

void
TpetraLinearSolver::setSystemObjects(
      Teuchos::RCP<LinSys::Matrix> matrix,
      Teuchos::RCP<LinSys::Vector> rhs)
{
  ThrowRequire(!matrix.is_null());
  ThrowRequire(!rhs.is_null());
  //ThrowRequire(solver_);

  matrix_ = matrix;
  rhs_ = rhs;
}

void TpetraLinearSolver::setupLinearSolver(
  Teuchos::RCP<LinSys::Vector> sln,
  Teuchos::RCP<LinSys::Matrix> matrix,
  Teuchos::RCP<LinSys::Vector> rhs,
  Teuchos::RCP<LinSys::MultiVector> coords)
{

  setSystemObjects(matrix,rhs);
  problem_ = Teuchos::RCP<LinSys::LinearProblem>(new LinSys::LinearProblem(matrix_, sln, rhs_) );

  if(activateMueLu_) {
    coords_ = coords;
  }
  else {
    Ifpack2::Factory factory;
    const std::string preconditionerType ("RELAXATION");
    preconditioner_ = factory.create (preconditionerType, Teuchos::rcp_const_cast<const LinSys::Matrix>(matrix_), 0);
    preconditioner_->setParameters(*paramsPrecond_);
    preconditioner_->initialize();
    problem_->setRightPrec(preconditioner_);

    // create the solver, e.g., gmres, cg, tfqmr, bicgstab
    LinSys::SolverFactory sFactory;
    solver_ = sFactory.create(config_->get_method(), params_);
    solver_->setProblem(problem_);
  }
}

void TpetraLinearSolver::destroyLinearSolver()
{
  problem_ = Teuchos::null;
  preconditioner_ = Teuchos::null;
  solver_ = Teuchos::null;
  coords_ = Teuchos::null;
  if (activateMueLu_) mueluPreconditioner_ = Teuchos::null;
}

void TpetraLinearSolver::setMueLu()
{
  if (solver_ != Teuchos::null && !recomputePreconditioner_ && !reusePreconditioner_) return;

  {
    Teuchos::RCP<Teuchos::Time> tm = Teuchos::TimeMonitor::getNewTimer("nalu MueLu preconditioner setup");
    Teuchos::TimeMonitor timeMon(*tm);

    if (recomputePreconditioner_ || mueluPreconditioner_ == Teuchos::null)
    {
      std::string xmlFileName = config_->muelu_xml_file();
      mueluPreconditioner_ = MueLu::CreateTpetraPreconditioner<SC,LO,GO,NO>(Teuchos::RCP<Tpetra::Operator<SC,LO,GO,NO> >(matrix_), xmlFileName, coords_);
    }
    else if (reusePreconditioner_) {
      MueLu::ReuseTpetraPreconditioner(matrix_, *mueluPreconditioner_);
    }
    if (config_->getSummarizeMueluTimer())
      Teuchos::TimeMonitor::summarize(std::cout, false, true, false, Teuchos::Union);
  }

  problem_->setRightPrec(mueluPreconditioner_);

  // create the solver, e.g., gmres, cg, tfqmr, bicgstab
  LinSys::SolverFactory sFactory;
  solver_ = sFactory.create(config_->get_method(), params_);
  solver_->setProblem(problem_);
}

int TpetraLinearSolver::residual_norm(int whichNorm, Teuchos::RCP<LinSys::Vector> sln, double& norm)
{
  LinSys::Vector resid(rhs_->getMap());
  ThrowRequire(! (sln.is_null()  || rhs_.is_null() ) );

  if (matrix_->isFillActive() )
  {
    // FIXME
    //!matrix_->fillComplete(map_, map_);
    throw std::runtime_error("residual_norm");
  }
  matrix_->apply(*sln, resid);

  LinSys::OneDVector rhs = rhs_->get1dViewNonConst ();
  LinSys::OneDVector res = resid.get1dViewNonConst ();
  for (int i=0; i<rhs.size(); ++i)
    res[i] -= rhs[i];
  if ( whichNorm == 0 )
    norm = resid.normInf();
  else if ( whichNorm == 1 )
    norm = resid.norm1();
  else if ( whichNorm == 2 )
    norm = resid.norm2();
  else
    return 1;

  return 0;
}


int
TpetraLinearSolver::solve(
  Teuchos::RCP<LinSys::Vector> sln,
  int & iters,
  double & finalResidNrm)
{
  ThrowRequire(!sln.is_null());


  const int status = 0;
  int whichNorm = 2;
  finalResidNrm=0.0;

  double time = -stk::cpu_time();
  if (activateMueLu_)
  {
    setMueLu();
  }
  else
  {
    preconditioner_->compute();
  }
  time += stk::cpu_time();
  timerPrecond_ += time;

  problem_->setProblem();
  solver_->solve();

  iters = solver_->getNumIters();
  residual_norm(whichNorm, sln, finalResidNrm);

  return status;
}

} // namespace nalu
} // namespace Sierra
