/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef LinearSolverConfig_h
#define LinearSolverConfig_h

#include <string>
#include <AztecOO.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

namespace YAML {
class Node;
}

namespace sierra{
namespace nalu{

class EpetraLinearSolverConfig {
  public:
    EpetraLinearSolverConfig();
    ~EpetraLinearSolverConfig() {}
    std::string name() const;
    void load(const YAML::Node & node);
    const int * aztec_options() const;
    const double * aztec_parameters() const;
    const Teuchos::RCP<Teuchos::ParameterList> ml_parameters() const;
    bool use_ml() const;
    bool use_mueLu() const { return useMueLu_; }
    std::string & muelu_xml_file() {return muelu_xml_file_;}
    bool getWriteMatrixFiles() { return writeMatrixFiles_; }
    bool getSummarizeMueluTimer() { return summarizeMueluTimer_; }
    bool recomputePreconditioner() { return recomputePreconditioner_; }
    bool reusePreconditioner() { return reusePreconditioner_; }
  private:
    static int string_to_AzSolver(const std::string & method);
    static int string_to_AzPrecond(const std::string & precond);
    static int string_to_AzSubdomainSolver(const std::string & solver);
  private:
    std::string name_;
    std::string method_;
    std::string precond_;
    std::string subMethod_;
    int az_options[AZ_OPTIONS_SIZE];
    double az_params[AZ_PARAMS_SIZE];
    bool writeMatrixFiles_;
    bool summarizeMueluTimer_;
    bool useML_;
    bool useMueLu_;
    Teuchos::RCP<Teuchos::ParameterList> mlParameterList_;
    std::string muelu_xml_file_;

    bool recomputePreconditioner_;
    bool reusePreconditioner_;
};

class TpetraLinearSolverConfig {
  public:
    TpetraLinearSolverConfig();
    ~TpetraLinearSolverConfig();
    std::string name() const;
    void load(const YAML::Node & node);
    const Teuchos::RCP<Teuchos::ParameterList> & params() const;
    const Teuchos::RCP<Teuchos::ParameterList> & paramsPrecond() const;
    bool getWriteMatrixFiles() { return writeMatrixFiles_; }
    bool getSummarizeMueluTimer() { return summarizeMueluTimer_; }
    bool use_MueLu() const {return useMueLu_;}
    std::string & muelu_xml_file() {return muelu_xml_file_;}
    bool recomputePreconditioner() { return recomputePreconditioner_; }
    bool reusePreconditioner() { return reusePreconditioner_; }
    std::string get_method() {return method_;}

  private:
    std::string name_;
    std::string method_;
    std::string precond_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
    Teuchos::RCP<Teuchos::ParameterList> paramsPrecond_;
    bool writeMatrixFiles_;
    bool summarizeMueluTimer_;

    std::string muelu_xml_file_;
    bool useMueLu_;

    bool recomputePreconditioner_;
    bool reusePreconditioner_;

};

} // namespace nalu
} // namespace Sierra

#endif
