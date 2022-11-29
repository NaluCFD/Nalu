/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef LinearSolverConfig_h
#define LinearSolverConfig_h

#include <string>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

namespace Ifpack2 {
class FunctionParameter;
}

namespace YAML {
class Node;
}

namespace sierra{
namespace nalu{

class LinearSolverConfig
{
public:
  LinearSolverConfig();
  virtual ~LinearSolverConfig() = default;

  virtual void load(const YAML::Node&) = 0;

  inline std::string name() const
  { return name_ ; }

  const Teuchos::RCP<Teuchos::ParameterList> & params() const
  { return params_; }

  const Teuchos::RCP<Teuchos::ParameterList> & paramsPrecond() const
  { return paramsPrecond_; }

  inline bool getWriteMatrixFiles() const
  { return writeMatrixFiles_; }

  inline bool recomputePreconditioner() const
  { return recomputePreconditioner_; }

  inline bool reusePreconditioner() const
  { return reusePreconditioner_; }

  std::string get_method() const
  {return method_;}

  std::string preconditioner_type() const
  { return preconditionerType_;}

  inline double tolerance() const { return tolerance_; }
  inline double finalTolerance() const { return finalTolerance_; }

  std::string solver_type() const
  { return solverType_; }

protected:
  std::string solverType_;
  std::string name_;
  std::string method_;
  std::string precond_;
  std::string preconditionerType_{"RELAXATION"};
  double tolerance_;
  double finalTolerance_;


  Teuchos::RCP<Teuchos::ParameterList> params_;
  Teuchos::RCP<Teuchos::ParameterList> paramsPrecond_;

  bool recomputePreconditioner_{true};
  bool reusePreconditioner_{false};
  bool writeMatrixFiles_{false};
};

class TpetraLinearSolverConfig : public LinearSolverConfig
{
public:
  TpetraLinearSolverConfig();
  virtual ~TpetraLinearSolverConfig();

  virtual void load(const YAML::Node & node) final;
  bool getSummarizeMueluTimer() { return summarizeMueluTimer_; }
  bool use_MueLu() const {return useMueLu_;}

private:
  bool summarizeMueluTimer_{false};
  bool useMueLu_{false};
};

} // namespace nalu
} // namespace Sierra

#endif
