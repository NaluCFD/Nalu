/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EpetraLinearSystem_h
#define EpetraLinearSystem_h

#include <LinearSystem.h>

#include <vector>
#include <string>

class Epetra_FECrsGraph;
class Epetra_FECrsMatrix;
class Epetra_Map;
class Epetra_Vector;
class Epetra_FEVector;

namespace sierra{
namespace nalu{

class Realm;
class LinearSolver;

class EpetraLinearSystem : public LinearSystem
{
public:

  EpetraLinearSystem(
    Realm &realm,
    const unsigned numDof,
    const std::string & name,
    LinearSolver * linearSolver);
  ~EpetraLinearSystem();

  // Graph/Matrix Construction
  void buildNodeGraph(const stk::mesh::PartVector & parts); // for nodal assembly (e.g., lumped mass and source)
  void buildFaceToNodeGraph(const stk::mesh::PartVector & parts); // face->node assembly
  void buildEdgeToNodeGraph(const stk::mesh::PartVector & parts); // edge->node assembly
  void buildElemToNodeGraph(const stk::mesh::PartVector & parts); // elem->node assembly
  void buildReducedElemToNodeGraph(const stk::mesh::PartVector & parts); // elem (nearest nodes only)->node assembly
  void buildFaceElemToNodeGraph(const stk::mesh::PartVector & parts); // elem:face->node assembly
  void buildEdgeHaloNodeGraph(const stk::mesh::PartVector & parts); // haloNode->elem_node assembly
  void buildNonConformalNodeGraph(const stk::mesh::PartVector & parts); // nonConformal->node assembly
  void buildOversetNodeGraph(const stk::mesh::PartVector & parts); // overset->elem_node assembly
  void finalizeLinearSystem();

  // Matrix Assembly
  void zeroSystem();

  void sumInto(
    const std::vector<stk::mesh::Entity> & sym_meshobj,
    std::vector<int> &scratchIds,
    std::vector<double> &scratchVals,
    const std::vector<double> & rhs,
    const std::vector<double> & lhs,
    const char *trace_tag=0
    );

  void applyDirichletBCs(
    stk::mesh::FieldBase * solutionField,
    stk::mesh::FieldBase * bcValuesField,
    const stk::mesh::PartVector & parts,
    const unsigned beginPos,
    const unsigned endPos);

  void prepareConstraints(
    const unsigned beginPos,
    const unsigned endPos) {}

  // Solve
  int solve(stk::mesh::FieldBase * linearSolutionField);
  void loadComplete();
  void writeToFile(const char * filename, bool useOwned=true);
  void writeSolutionToFile(const char * filename, bool useOwned=true);

private:
  void beginLinearSystemConstruction();
  void checkError(
    const int err_code,
    const char * msg);

  void copy_epetra_to_stk(
    const Epetra_Vector * epetraField,
    stk::mesh::FieldBase * stkField);

  Epetra_Map *rowMap_;
  Epetra_FECrsGraph * graph_;
  Epetra_FECrsMatrix * lhs_;
  Teuchos::RCP<Epetra_FECrsMatrix> mueLuLHS;
  Epetra_FEVector * rhs_;
  Epetra_Vector * sln_;
  double * coords_;

};


} // namespace nalu
} // namespace Sierra

#endif
