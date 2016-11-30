/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <LinearSolverConfig.h>
#include <NaluEnv.h>
#include <NaluParsing.h>
#include <yaml-cpp/yaml.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <ml_MultiLevelPreconditioner.h>
#include <BelosTypes.hpp>

#include <ostream>

namespace sierra{
namespace nalu{

TpetraLinearSolverConfig::TpetraLinearSolverConfig() :
  params_(Teuchos::rcp(new Teuchos::ParameterList)),
  paramsPrecond_(Teuchos::rcp(new Teuchos::ParameterList)),
  useMueLu_(false),
  recomputePreconditioner_(true),
  reusePreconditioner_(false),
  writeMatrixFiles_(false),
  summarizeMueluTimer_(false),
  preconditionerType_("RELAXATION")
{}

TpetraLinearSolverConfig::~TpetraLinearSolverConfig()
{
}

std::string
TpetraLinearSolverConfig::name() const
{
  return name_;
}

const Teuchos::RCP<Teuchos::ParameterList> &
TpetraLinearSolverConfig::params() const
{
  return params_;
}

const Teuchos::RCP<Teuchos::ParameterList> &
TpetraLinearSolverConfig::paramsPrecond() const
{
  return paramsPrecond_;
}

void
TpetraLinearSolverConfig::load(const YAML::Node & node)
{
  name_ = node["name"].as<std::string>() ;
  method_ = node["method"].as<std::string>() ;
  get_if_present(node, "preconditioner", precond_, std::string("default"));

  double tol;
  int max_iterations, kspace, output_level;

  get_if_present(node, "tolerance", tol, 1.e-4);
  get_if_present(node, "max_iterations", max_iterations, 50);
  get_if_present(node, "kspace", kspace, 50);
  get_if_present(node, "output_level", output_level, 0);

  //Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::params();
  params_->set("Convergence Tolerance", tol);
  params_->set("Maximum Iterations", max_iterations);
  if (output_level > 0)
  {
    params_->set("Verbosity", Belos::Debug + Belos::Warnings + Belos::IterationDetails
      + Belos::OrthoDetails + Belos::FinalSummary
      + Belos::TimingDetails + Belos::StatusTestDetails);
  }

  params_->set("Output Frequency", output_level);
  Teuchos::RCP<std::ostream> belosOutputStream = Teuchos::rcpFromRef (NaluEnv::self().naluOutputP0());
  params_->set("Output Stream", belosOutputStream);
  params_->set("Num Blocks", kspace);
  params_->set("Maximum Restarts", std::max(1,max_iterations/kspace));
  std::string orthoType = "ICGS";
  params_->set("Orthogonalization",orthoType);
  params_->set("Implicit Residual Scaling", "Norm of Preconditioned Initial Residual");

  if (precond_ == "sgs") {
    preconditionerType_ = "RELAXATION";
    paramsPrecond_->set("relaxation: type","Symmetric Gauss-Seidel");
    paramsPrecond_->set("relaxation: sweeps",1);
  }
  else if (precond_ == "jacobi" || precond_ == "default") {
    preconditionerType_ = "RELAXATION";
    paramsPrecond_->set("relaxation: type","Jacobi");
    paramsPrecond_->set("relaxation: sweeps",1);
  }
  else if (precond_ == "ilut" ) {
    preconditionerType_ = "ILUT";
  }
  else if (precond_ == "riluk" ) {
    preconditionerType_ = "RILUK";
  }
  else if (precond_ == "muelu") {
    muelu_xml_file_ = std::string("milestone.xml");
    get_if_present(node, "muelu_xml_file_name", muelu_xml_file_, muelu_xml_file_);
    useMueLu_ = true;
  }
  else {
    throw std::runtime_error("invalid linear solver preconditioner specified ");
  }

  get_if_present(node, "write_matrix_files", writeMatrixFiles_, writeMatrixFiles_);
  get_if_present(node, "summarize_muelu_timer", summarizeMueluTimer_, summarizeMueluTimer_);

  get_if_present(node, "recompute_preconditioner", recomputePreconditioner_, recomputePreconditioner_);
  get_if_present(node, "reuse_preconditioner",     reusePreconditioner_,     reusePreconditioner_);

}

} // namespace nalu
} // namespace Sierra
