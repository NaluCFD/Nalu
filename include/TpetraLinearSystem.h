/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TpetraLinearSystem_h
#define TpetraLinearSystem_h

#include <LinearSystem.h>

#include <Tpetra_DefaultPlatform.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <stk_mesh/base/Entity.hpp>

#include <vector>
#include <string>
#include <boost/unordered_map.hpp>

namespace stk {
namespace mesh {
  typedef uint64_t EntityId;
}
}

namespace sierra{
namespace nalu{

class Realm;
class LinearSolver;

typedef boost::unordered_map<stk::mesh::EntityId, size_t>  MyLIDMapType;

class TpetraLinearSystem : public LinearSystem
{
public:
  typedef LinSys::GlobalOrdinal GlobalOrdinal;
  typedef LinSys::LocalOrdinal  LocalOrdinal;

  TpetraLinearSystem(
    Realm &realm,
    const unsigned numDof,
    const std::string & name,
    LinearSolver * linearSolver);
  ~TpetraLinearSystem();

   // Graph/Matrix Construction
  void buildNodeGraph(const stk::mesh::PartVector & parts); // for nodal assembly (e.g., lumped mass and source)
  void buildFaceToNodeGraph(const stk::mesh::PartVector & parts); // face->node assembly
  void buildEdgeToNodeGraph(const stk::mesh::PartVector & parts); // edge->node assembly
  void buildElemToNodeGraph(const stk::mesh::PartVector & parts); // elem->node assembly
  void buildReducedElemToNodeGraph(const stk::mesh::PartVector & parts); // elem (nearest nodes only)->node assembly
  void buildFaceElemToNodeGraph(const stk::mesh::PartVector & parts); // elem:face->node assembly
  void buildNonConformalNodeGraph(const stk::mesh::PartVector & parts); // nonConformal->node assembly
  void buildOversetNodeGraph(const stk::mesh::PartVector & parts); // overset->elem_node assembly
  void finalizeLinearSystem();

  // Matrix Assembly
  void zeroSystem();

  void sumInto(
    const std::vector<stk::mesh::Entity> & entities,
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
    const unsigned endPos);

  // Solve
  int solve(stk::mesh::FieldBase * linearSolutionField);
  void loadComplete();
  void writeToFile(const char * filename, bool useOwned=true);
  void printInfo(bool useOwned=true);
  void writeSolutionToFile(const char * filename, bool useOwned=true);
  size_t lookup_myLID(MyLIDMapType& myLIDs, stk::mesh::EntityId entityId, const char* msg=nullptr, stk::mesh::Entity entity = stk::mesh::Entity());

  enum DOFStatus {
    DS_NotSet           = 0,
    DS_SkippedDOF       = 1 << 1,
    DS_OwnedDOF         = 1 << 2,
    DS_GloballyOwnedDOF = 1 << 3,
    DS_GhostedDOF       = 1 << 4
  };

  int getDofStatus(stk::mesh::Entity node);

private:
  void beginLinearSystemConstruction();

  void checkError(
    const int err_code,
    const char * msg) {}

  void copy_tpetra_to_stk(
    const Teuchos::RCP<LinSys::Vector> tpetraVector,
    stk::mesh::FieldBase * stkField);

  // This method copies a stk::mesh::field to a tpetra multivector. Each dof/node is written into a different
  // vector in the multivector.
  void copy_stk_to_tpetra(stk::mesh::FieldBase * stkField,
    const Teuchos::RCP<LinSys::MultiVector> tpetraVector);

  void addConnections(const stk::mesh::Entity* entities, size_t num_entities);
  void checkForNaN(bool useOwned);
  bool checkForZeroRow(bool useOwned, bool doThrow, bool doPrint=false);

  typedef std::pair<stk::mesh::Entity, stk::mesh::Entity> Connection;
  typedef std::set< Connection > ConnectionSet;
  typedef std::vector< Connection > ConnectionVec;
  ConnectionSet connectionSet_;
  std::vector<GlobalOrdinal> totalGids_;

  Teuchos::RCP<LinSys::Node>   node_;

  // all rows, otherwise known as col map
  Teuchos::RCP<LinSys::Map>    totalColsMap_;

  // Map of rows my proc owns (locally owned)
  Teuchos::RCP<LinSys::Map>    ownedRowsMap_;

  // Map of all rows my proc references
  Teuchos::RCP<LinSys::Map>    ownedPlusGloballyOwnedRowsMap_;

  // Only nodes that share with other procs that I don't own = Global = !locally owned
  Teuchos::RCP<LinSys::Map>    globallyOwnedRowsMap_;

  Teuchos::RCP<LinSys::Graph>  ownedGraph_;
  Teuchos::RCP<LinSys::Graph>  globallyOwnedGraph_;

  Teuchos::RCP<LinSys::Matrix> ownedMatrix_;
  Teuchos::RCP<LinSys::Vector> ownedRhs_;

  Teuchos::RCP<LinSys::Matrix> globallyOwnedMatrix_;
  Teuchos::RCP<LinSys::Vector> globallyOwnedRhs_;

  Teuchos::RCP<LinSys::Vector> sln_;
  Teuchos::RCP<LinSys::Vector> globalSln_;
  Teuchos::RCP<LinSys::Export> exporter_;
  Teuchos::RCP<LinSys::Import> importer_;

  MyLIDMapType myLIDs_;
  LocalOrdinal maxOwnedRowId_; // = num_owned_nodes * numDof_
  LocalOrdinal maxGloballyOwnedRowId_; // = (num_owned_nodes + num_globallyOwned_nodes) * numDof_
};


} // namespace nalu
} // namespace Sierra

#endif
