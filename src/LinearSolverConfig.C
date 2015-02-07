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

EpetraLinearSolverConfig::EpetraLinearSolverConfig() :
  useML_(false),
  useMueLu_(false),
  mlParameterList_(Teuchos::rcp(new Teuchos::ParameterList))
{
  for (int i = 0 ; i < AZ_OPTIONS_SIZE; ++i ) az_options[i] = 0;
  for (int i = 0 ; i < AZ_PARAMS_SIZE; ++i )  az_params[i] = 0.0;
}

std::string
EpetraLinearSolverConfig::name() const
{
  return name_;
}

const int *
EpetraLinearSolverConfig::aztec_options() const
{
  return az_options;
}

const double *
EpetraLinearSolverConfig::aztec_parameters() const
{
  return az_params;
}

const Teuchos::RCP<Teuchos::ParameterList>
EpetraLinearSolverConfig::ml_parameters() const
{
  return mlParameterList_;
}

bool EpetraLinearSolverConfig::use_ml() const {return useML_;}

void
EpetraLinearSolverConfig::load(const YAML::Node & node)
{
  AZ_defaults(az_options, az_params);
  node["name"]   >> name_;
  node["method"] >> method_;
  get_if_present(node, "preconditioner", precond_, std::string("default"));
  get_if_present(node, "subdomain_solver", subMethod_, std::string("default"));
  if (precond_ == "ML")
  {
    useML_ = true;

    ML_Epetra::SetDefaults("SA", *mlParameterList_);

    const YAML::Node * int_nodes = node.FindValue("ML_options_int");
    if ( int_nodes )
    {
      for ( size_t inode = 0; inode <  int_nodes->size(); ++inode )
      {
        const YAML::Node & integer_parameter_node = (* int_nodes)[inode];
        std::string option_name;
        int option_value;
        integer_parameter_node["name"] >> option_name;
        integer_parameter_node["value"] >> option_value;
        mlParameterList_->set(option_name,option_value);

      }
    }

    const YAML::Node * str_nodes = node.FindValue("ML_options_string");
    if ( str_nodes )
    {
      for ( size_t inode = 0; inode <  str_nodes->size(); ++inode )
      {
        const YAML::Node & integer_parameter_node = (* str_nodes)[inode];
        std::string option_name;
        std::string option_value;
        integer_parameter_node["name"] >> option_name;
        integer_parameter_node["value"] >> option_value;
        mlParameterList_->set(option_name,option_value);

      }
    }

    const YAML::Node * real_nodes = node.FindValue("ML_options_real");
    if ( real_nodes )
    {
      for ( size_t inode = 0; inode <  real_nodes->size(); ++inode )
      {
        const YAML::Node & integer_parameter_node = (* real_nodes)[inode];
        std::string option_name;
        double option_value;
        integer_parameter_node["name"] >> option_name;
        integer_parameter_node["value"] >> option_value;
        mlParameterList_->set(option_name,option_value);

      }
    }

  }
  else if (precond_ == "muelu")
  {
    useMueLu_ = true;
    muelu_xml_file_ = std::string("milestone.xml");
    get_if_present(node, "muelu_xml_file_name", muelu_xml_file_, muelu_xml_file_);
    az_options[AZ_precond] = string_to_AzPrecond("multilevel");
  }
  else
  {
    az_options[AZ_precond] = string_to_AzPrecond(precond_);
  }
  az_options[AZ_solver] = string_to_AzSolver(method_);
  az_options[AZ_subdomain_solve] = string_to_AzSubdomainSolver(subMethod_);
  az_options[AZ_conv] = AZ_r0; // should be the default
  get_if_present(node, "output_level", az_options[AZ_output], 0);

  get_if_present_no_default(node, "max_iterations", az_options[AZ_max_iter]);
  get_if_present_no_default(node, "kspace", az_options[AZ_kspace]);
  get_if_present_no_default(node, "graph_fill", az_options[AZ_graph_fill]);
  get_if_present_no_default(node, "residual_norm", az_options[AZ_resid]);
  get_if_present_no_default(node, "overlap", az_options[AZ_overlap]);
  get_if_present_no_default(node, "type_overlap", az_options[AZ_type_overlap]);
  get_if_present_no_default(node, "poly_ord", az_options[AZ_poly_ord]);

  get_if_present_no_default(node, "tolerance", az_params[AZ_tol]);
  get_if_present_no_default(node, "drop_tolerance", az_params[AZ_drop]);
  get_if_present_no_default(node, "ilut_fill", az_params[AZ_ilut_fill]);
  get_if_present_no_default(node, "omega", az_params[AZ_omega]);
  
  get_if_present(node, "write_matrix_files", writeMatrixFiles_, false);
  get_if_present(node, "summarize_muelu_timer", summarizeMueluTimer_, false);
  
  get_if_present(node, "recompute_preconditioner", recomputePreconditioner_, true);
  get_if_present(node, "reuse_preconditioner",     reusePreconditioner_,     false);
}

int
EpetraLinearSolverConfig::string_to_AzSolver(const std::string & method)
{
  if(method == "default")  return AZ_default;
  if(method == "cg")       return AZ_cg;
  if(method == "gmres")    return AZ_gmres;
  if(method == "cgs")      return AZ_cgs;
  if(method == "tfqmr")    return AZ_tfqmr;
  if(method == "bicgstab") return AZ_bicgstab;
  if(method == "slu")      return AZ_slu;
  if(method == "symmlq")   return AZ_symmlq;
  if(method == "gmresr")   return AZ_GMRESR;
  if(method == "fixed_pt") return AZ_fixed_pt;
  if(method == "analyze")  return AZ_analyze;
  if(method == "lu")       return AZ_lu;
  throw std::runtime_error("invalid linear solver method specified ");
}

int
EpetraLinearSolverConfig::string_to_AzPrecond(const std::string & precond)
{
  if(precond == "default")      return AZ_default;
  if(precond == "none")         return AZ_none;
  if(precond == "jacobi")       return AZ_Jacobi;
  if(precond == "sgs")          return AZ_sym_GS;
  if(precond == "neumann")      return AZ_Neumann;
  if(precond == "ls")           return AZ_ls;
  if(precond == "ilu")          return AZ_ilu;
  if(precond == "bilu")         return AZ_bilu;
  if(precond == "icc")          return AZ_icc;
  if(precond == "ilut")         return AZ_ilut;
  if(precond == "rilu")         return AZ_rilu;
  if(precond == "recursive")    return AZ_recursive;
  if(precond == "smoother")     return AZ_smoother;
  if(precond == "dom_decomp")   return AZ_dom_decomp;
  if(precond == "multilevel")   return AZ_multilevel;
  if(precond == "user_solve")   return AZ_user_precond;
  if(precond == "bilu_ifp")     return AZ_bilu_ifp;
  throw std::runtime_error("invalid linear solver preconditioner specified ");
}

int
EpetraLinearSolverConfig::string_to_AzSubdomainSolver(const std::string & solver)
{
  if(solver == "default")      return AZ_default;
  if(solver == "none")         return AZ_none;
  if(solver == "jacobi")       return AZ_Jacobi;
  if(solver == "sgs")          return AZ_sym_GS;
  if(solver == "neumann")      return AZ_Neumann;
  if(solver == "ls")           return AZ_ls;
  if(solver == "ilu")          return AZ_ilu;
  if(solver == "bilu")         return AZ_bilu;
  if(solver == "icc")          return AZ_icc;
  if(solver == "ilut")         return AZ_ilut;
  if(solver == "rilu")         return AZ_rilu;
  if(solver == "recursive")    return AZ_recursive;
  if(solver == "smoother")     return AZ_smoother;
  if(solver == "dom_decomp")   return AZ_dom_decomp;
  if(solver == "multilevel")   return AZ_multilevel;
  if(solver == "user_solve")   return AZ_user_precond;
  if(solver == "bilu_ifp")     return AZ_bilu_ifp;
  throw std::runtime_error("invalid linear subdomain solver specified ");
}

TpetraLinearSolverConfig::TpetraLinearSolverConfig() :
  params_(Teuchos::rcp(new Teuchos::ParameterList)),
  paramsPrecond_(Teuchos::rcp(new Teuchos::ParameterList)),
  useMueLu_(false)
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
  node["name"]    >> name_;
  node["method"]  >> method_;
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
    paramsPrecond_->set("relaxation: type","Symmetric Gauss-Seidel");
    paramsPrecond_->set("relaxation: sweeps",1);
  }
  else if (precond_ == "jacobi" || precond_ == "default") {
    paramsPrecond_->set("relaxation: type","Jacobi");
    paramsPrecond_->set("relaxation: sweeps",1);
  }
  else if (precond_ == "muelu") {
    muelu_xml_file_ = std::string("milestone.xml");
    get_if_present(node, "muelu_xml_file_name", muelu_xml_file_, muelu_xml_file_);
    useMueLu_ = true;
  }
  else {
    throw std::runtime_error("invalid linear solver preconditioner specified ");
  }

  get_if_present(node, "write_matrix_files", writeMatrixFiles_, false);
  get_if_present(node, "summarize_muelu_timer", summarizeMueluTimer_, false);

  get_if_present(node, "recompute_preconditioner", recomputePreconditioner_, true);
  get_if_present(node, "reuse_preconditioner",     reusePreconditioner_,     false);

}

} // namespace nalu
} // namespace Sierra
