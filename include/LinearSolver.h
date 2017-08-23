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

/** Type of solvers available in Nalu simulation **/
  enum PetraType {
    PT_TPETRA,       //!< Nalu Tpetra interface
    PT_HYPRE,        //!< Direct HYPRE interface
    PT_TPETRA_HYPRE, //!< Tpetra to Hypre interface via xSDK
    PT_END
  };


class LinearSolvers;
class Simulation;
class Realm;

/** An abstract representation of a linear solver in Nalu
 *
 *  Defines the basic API supported by the linear solvers for use within Nalu.
 *  See concrete implementations such as sierra::nalu::TpetraLinearSolver for
 *  more details.
 */
class LinearSolver
{
  public:
    LinearSolver(
      std::string name,
      LinearSolvers* linearSolvers,
      LinearSolverConfig* config)
      : name_(name),
        linearSolvers_(linearSolvers),
        config_(config),
        recomputePreconditioner_(config->recomputePreconditioner()),
        reusePreconditioner_(config->reusePreconditioner()),
        timerPrecond_(0.0)
    {}
    virtual ~LinearSolver() {}

  //! User-friendly identifier for this particular solver instance
    std::string name_;

  //! Type of solver instance as defined in sierra::nalu::PetraType
    virtual PetraType getType() = 0;

  /** Solve the linear system Ax = b
   *
   *  @param[out] sln The solution vector
   *  @param[out] iters The number of linear solver iterations to convergence
   *  @param[out] finalResidNrm The final residual norm
   *  @param[in]  isFinalOuterIter Is this the final outer iteration
   */
    virtual int solve(Teuchos::RCP<LinSys::Vector>, int&, double&, bool) = 0;

  /** Create the linear solver instances for Trilinos Belos solvers
   *
   * @param[in] sln The solution vector instance
   * @param[in] matrix The Tpetra matrix instance
   * @param[in] rhs The Tpetra vector instance for RHS
   * @param[in] coords The grid coordinates as a Tpetra data structure
   */
    virtual void setupLinearSolver(
      Teuchos::RCP<LinSys::Vector>,
      Teuchos::RCP<LinSys::Matrix>,
      Teuchos::RCP<LinSys::Vector>,
      Teuchos::RCP<LinSys::MultiVector>) = 0;

  /** Utility method to cleanup solvers during simulation
   */
    virtual void destroyLinearSolver() = 0;

    Simulation* root();
    LinearSolvers* parent();
    LinearSolvers* linearSolvers_;
    Realm* realm_{nullptr};

  //! The number of degrees of freedom in the equation system. Default: 1
    int numDof_{1};

  protected:
  LinearSolverConfig* config_;
  bool recomputePreconditioner_;
  bool reusePreconditioner_;
  double timerPrecond_;
  bool activateMueLu_{false};

  public:
  //! Flag indicating whether the preconditioner is recomputed on each invocation
  bool & recomputePreconditioner() {return recomputePreconditioner_;}
  //! Flag indicating whether the preconditioner is reused on each invocation
  bool & reusePreconditioner() {return reusePreconditioner_;}

  //! Reset the preconditioner timer to 0.0 for future accumulation
  void zero_timer_precond() { timerPrecond_ = 0.0;}

  //! Get the preconditioner timer for the last invocation
  double get_timer_precond() { return timerPrecond_;}

  //! Flag indicating whether the user has activated MueLU
  bool& activeMueLu() { return activateMueLu_; }

  //! Get the solver configuration specified in the input file
  LinearSolverConfig* getConfig() { return config_; }
};

class TpetraLinearSolver : public LinearSolver
{
  public:

  /**
   *  @param[in] solverName The name of the solver
   *  @param[in] config Solver configuration
   */
  TpetraLinearSolver(
    std::string solverName,
    TpetraLinearSolverConfig *config,
    const Teuchos::RCP<Teuchos::ParameterList> params,
    const Teuchos::RCP<Teuchos::ParameterList> paramsPrecond,
    LinearSolvers *linearSolvers);
  virtual ~TpetraLinearSolver() ;
  
    void setSystemObjects(
      Teuchos::RCP<LinSys::Matrix> matrix,
      Teuchos::RCP<LinSys::Vector> rhs);

    virtual void setupLinearSolver(
      Teuchos::RCP<LinSys::Vector> sln,
      Teuchos::RCP<LinSys::Matrix> matrix,
      Teuchos::RCP<LinSys::Vector> rhs,
      Teuchos::RCP<LinSys::MultiVector> coords) override;

    virtual void destroyLinearSolver() override;

  //! Initialize the MueLU preconditioner before solve
    void setMueLu();

  /** Compute the norm of the non-linear solution vector
   *
   *  @param[in] whichNorm [0, 1, 2] norm to be computed
   *  @param[in] sln The solution vector
   *  @param[out] norm The norm of the solution vector
   */
    int residual_norm(int whichNorm, Teuchos::RCP<LinSys::Vector> sln, double& norm);

  /** Solve the linear system Ax = b
   *
   *  @param[out] sln The solution vector
   *  @param[out] iters The number of linear solver iterations to convergence
   *  @param[out] finalResidNrm The final residual norm
   *  @param[in]  isFinalOuterIter Is this the final outer iteration
   */
    virtual int solve(
      Teuchos::RCP<LinSys::Vector> sln,
      int & iterationCount,
      double & scaledResidual,
      bool isFinalOuterIter) override;

    virtual PetraType getType() override { return PT_TPETRA; }

  private:
  //! The solver parameters
    const Teuchos::RCP<Teuchos::ParameterList> params_;

  //! The preconditioner parameters
    const Teuchos::RCP<Teuchos::ParameterList> paramsPrecond_;
    Teuchos::RCP<LinSys::Matrix> matrix_;
    Teuchos::RCP<LinSys::Vector> rhs_;
    Teuchos::RCP<LinSys::LinearProblem> problem_;
    Teuchos::RCP<LinSys::SolverManager> solver_;
    Teuchos::RCP<LinSys::Preconditioner> preconditioner_;
    Teuchos::RCP<MueLu::TpetraOperator<SC,LO,GO,NO> > mueluPreconditioner_;
    Teuchos::RCP<LinSys::MultiVector> coords_;

    std::string preconditionerType_;
};

} // namespace nalu
} // namespace Sierra

#endif
