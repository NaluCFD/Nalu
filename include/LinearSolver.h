/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef LinearSolver_h
#define LinearSolver_h

#include <LinearSolverTypes.h>
#include <LinearSolverConfig.h>
#include <ml_MultiLevelPreconditioner.h>

#include <LinearSolverTypes.h>

#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#include <Ifpack2_Factory.hpp>

// Header files defining default types for template parameters.
// These headers must be included after other MueLu/Xpetra headers.
typedef double                                                        Scalar;
typedef long                                                          GlobalOrdinal;
typedef int                                                           LocalOrdinal;
typedef Tpetra::DefaultPlatform::DefaultPlatformType                  Platform;
typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal>::node_type           Node;
typedef Teuchos::ScalarTraits<Scalar> STS;

// MueLu main header: include most common header files in one line
#include <MueLu.hpp>

#include <MueLu_TrilinosSmoother.hpp> //TODO: remove
#include <MueLu_TpetraOperator.hpp>

#include <MueLu_UseShortNames.hpp>    // => typedef MueLu::FooClass<Scalar, LocalOrdinal, ...> Foo

namespace sierra{
namespace nalu{

  enum PetraType {
    PT_TPETRA,
    PT_END
  };


class LinearSolvers;
class Simulation;

class LinearSolver
{
  public:
  LinearSolver(std::string name, LinearSolvers *linearSolvers,
    bool recompute_preconditioner, bool reuse_preconditioner) : name_(name), linearSolvers_(linearSolvers),
    recomputePreconditioner_(recompute_preconditioner), reusePreconditioner_(reuse_preconditioner), timerPrecond_(0.0) {}
  virtual ~LinearSolver() {}
  std::string name_;
  virtual PetraType getType() = 0;
  Simulation *root();
  LinearSolvers *parent();
  LinearSolvers *linearSolvers_;
  protected:
  bool recomputePreconditioner_;
  bool reusePreconditioner_;
  double timerPrecond_;
  public:
  bool & recomputePreconditioner() {return recomputePreconditioner_;}
  bool & reusePreconditioner() {return reusePreconditioner_;}
  void zero_timer_precond() { timerPrecond_ = 0.0;}
  double get_timer_precond() { return timerPrecond_;}
};

class TpetraLinearSolver : public LinearSolver
{
  public:

  TpetraLinearSolver(
    std::string solverName,
    TpetraLinearSolverConfig *config,
    const Teuchos::RCP<Teuchos::ParameterList> params, const Teuchos::RCP<Teuchos::ParameterList> paramsPrecond,
    LinearSolvers *linearSolvers);
  ~TpetraLinearSolver() ;
  
    void setSystemObjects(
      Teuchos::RCP<LinSys::Matrix> matrix,
      Teuchos::RCP<LinSys::Vector> rhs);

    void setupLinearSolver(
      Teuchos::RCP<LinSys::Vector> sln,
      Teuchos::RCP<LinSys::Matrix> matrix,
      Teuchos::RCP<LinSys::Vector> rhs,
      Teuchos::RCP<LinSys::MultiVector> coords);

    void destroyLinearSolver();

    void setMueLu();

    int residual_norm(int whichNorm, Teuchos::RCP<LinSys::Vector> sln, double& norm);

    int solve(
      Teuchos::RCP<LinSys::Vector> sln,
      int & iterationCount,
      double & scaledResidual);

    virtual PetraType getType() { return PT_TPETRA; }
    TpetraLinearSolverConfig *getConfig() { return config_; }

    bool & activeMueLu(){ return activateMueLu_; }

  private:
    TpetraLinearSolverConfig *config_;
    const Teuchos::RCP<Teuchos::ParameterList> params_;
    const Teuchos::RCP<Teuchos::ParameterList> paramsPrecond_;
    Teuchos::RCP<LinSys::Matrix> matrix_;
    Teuchos::RCP<LinSys::Vector> rhs_;
    Teuchos::RCP<LinSys::LinearProblem> problem_;
    Teuchos::RCP<LinSys::SolverManager> solver_;
    Teuchos::RCP<LinSys::Preconditioner> preconditioner_;
    Teuchos::RCP<MueLu::TpetraOperator<SC,LO,GO,NO> > mueluPreconditioner_;
    Teuchos::RCP<LinSys::MultiVector> coords_;

    bool activateMueLu_;
    std::string preconditionerType_;
};

} // namespace nalu
} // namespace Sierra

#endif
