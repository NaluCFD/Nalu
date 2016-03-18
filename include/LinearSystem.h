/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef LinearSystem_h
#define LinearSystem_h

#include <LinearSolverTypes.h>

#include <Teuchos_RCP.hpp>
#include <Tpetra_DefaultPlatform.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#include <vector>
#include <string>

namespace stk { namespace mesh { struct Entity; } }

namespace stk{
namespace mesh{
class FieldBase;
class Part;
typedef std::vector< Part * > PartVector ;
}
}
namespace sierra{
namespace nalu{

class Realm;
class LinearSolver;

class LinearSystem
{
public:

  LinearSystem(
    Realm &realm,
    const unsigned numDof,
    const std::string & name,
    LinearSolver *linearSolver);

  virtual ~LinearSystem() {}

  static LinearSystem *create(Realm& realm, const unsigned numDof, const std::string & name, LinearSolver *linearSolver);

  // Graph/Matrix Construction
  virtual void buildNodeGraph(const stk::mesh::PartVector & parts)=0; // for nodal assembly (e.g., lumped mass and source)
  virtual void buildFaceToNodeGraph(const stk::mesh::PartVector & parts)=0; // face->node assembly
  virtual void buildEdgeToNodeGraph(const stk::mesh::PartVector & parts)=0; // edge->node assembly
  virtual void buildElemToNodeGraph(const stk::mesh::PartVector & parts)=0; // elem->node assembly
  virtual void buildReducedElemToNodeGraph(const stk::mesh::PartVector & parts)=0; // elem (nearest nodes only)->node assembly
  virtual void buildFaceElemToNodeGraph(const stk::mesh::PartVector & parts)=0; // elem:face->node assembly
  virtual void buildEdgeHaloNodeGraph(const stk::mesh::PartVector & parts)=0; // haloNode->elem_node assembly
  virtual void buildNonConformalNodeGraph(const stk::mesh::PartVector & parts)=0; // nonConformal->elem_node assembly
  virtual void buildOversetNodeGraph(const stk::mesh::PartVector & parts)=0; // overset->elem_node assembly
  virtual void finalizeLinearSystem()=0;

  // Matrix Assembly
  virtual void zeroSystem()=0;

  virtual void sumInto(
    const std::vector<stk::mesh::Entity> & sym_meshobj,
    std::vector<int> &scratchIds,
    std::vector<double> &scratchVals,
    const std::vector<double> & rhs,
    const std::vector<double> & lhs,
    const char *trace_tag=0
    )=0;

  virtual void applyDirichletBCs(
    stk::mesh::FieldBase * solutionField,
    stk::mesh::FieldBase * bcValuesField,
    const stk::mesh::PartVector & parts,
    const unsigned beginPos,
    const unsigned endPos)=0;

  virtual void prepareConstraints(
    const unsigned beginPos,
    const unsigned endPos)=0;

  // Solve
  virtual int solve(stk::mesh::FieldBase * linearSolutionField)=0;
  virtual void loadComplete()=0;

  virtual void writeToFile(const char * filename, bool useOwned=true)=0;
  virtual void writeSolutionToFile(const char * filename, bool useOwned=true)=0;
  unsigned numDof() const { return numDof_; }
  const int & linearSolveIterations() {return linearSolveIterations_; }
  const double & linearResidual() {return linearResidual_; }
  const double & nonLinearResidual() {return nonLinearResidual_; }
  const double & scaledNonLinearResidual() {return scaledNonLinearResidual_; }
  void setNonLinearResidual(const double nlr) { nonLinearResidual_ = nlr;}
  const std::string name() { return name_; }
  bool & recomputePreconditioner() {return recomputePreconditioner_;}
  bool & reusePreconditioner() {return reusePreconditioner_;}
  double get_timer_precond();
  void zero_timer_precond();

protected:
  virtual void beginLinearSystemConstruction()=0;
  virtual void checkError(
    const int err_code,
    const char * msg)=0;

  void sync_field(const stk::mesh::FieldBase *field);
  bool debug();

  Realm &realm_;
  bool inConstruction_;
  int writeCounter_;

  const unsigned numDof_;
  const std::string name_;
  LinearSolver * linearSolver_;
  int linearSolveIterations_;
  double nonLinearResidual_;
  double linearResidual_;
  double firstNonLinearResidual_;
  double scaledNonLinearResidual_;
  bool recomputePreconditioner_;
  bool reusePreconditioner_;

public:
  bool provideOutput_;
};

} // namespace nalu
} // namespace Sierra

#endif
