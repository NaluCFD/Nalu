/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef HYPRELINEARSOLVER_H
#define HYPRELINEARSOLVER_H

#include "LinearSolver.h"
#include "NaluHypre.h"

namespace sierra {
namespace nalu {

class HypreLinearSolver: public LinearSolver
{
public:
  HypreLinearSolver(
    std::string name,
    HypreLinearSolverConfig* config,
    LinearSolvers* linearSolvers);

  virtual ~HypreLinearSolver();

  virtual void setupLinearSolver(
    Teuchos::RCP<LinSys::Vector> sln,
    Teuchos::RCP<LinSys::Matrix> matrix,
    Teuchos::RCP<LinSys::Vector> rhs,
    Teuchos::RCP<LinSys::MultiVector> coords) override;

  virtual void destroyLinearSolver() override;

  virtual  int solve(
    Teuchos::RCP<LinSys::Vector> sln,
    int & iterationCount,
    double & scaledResidual,
    bool isFinalOuterIter) override;

  virtual PetraType getType() override { return PT_TPETRA_HYPRE; }

protected:
  void sync_hypre_id();

  Teuchos::RCP<LinSys::Matrix> matrix_;
  Teuchos::RCP<LinSys::Vector> rhs_;
  Teuchos::RCP<Ifpack2::NaluHypre<LinSys::Scalar, LinSys::LocalOrdinal, LinSys::GlobalOrdinal, LinSys::Node>> hypreSolver_;
};

}  // nalu
}  // sierra


#endif /* HYPRELINEARSOLVER_H */
