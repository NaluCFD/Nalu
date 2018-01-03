/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TpetraLinearSystem_h
#define TpetraLinearSystem_h

#include <LinearSystem.h>

#include <KokkosInterface.h>

#include <Tpetra_DefaultPlatform.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/FieldBase.hpp>

#include <vector>
#include <string>
#include <boost/unordered_map.hpp>

#include <Kokkos_UnorderedMap.hpp>

namespace sierra {
namespace nalu {

class Realm;
class EquationSystem;
class LinearSolver;

typedef boost::unordered_map<stk::mesh::EntityId, size_t>  MyLIDMapType;

typedef std::pair<stk::mesh::Entity, stk::mesh::Entity> Connection;
typedef Kokkos::UnorderedMap<Connection,void> ConnectionSetKK;
typedef std::vector< Connection > ConnectionVec;

typedef typename LinSys::Vector::dual_view_type dual_view_type;
typedef typename dual_view_type::t_host host_view_type;

  enum DOFStatus {
    DS_NotSet           = 0,
    DS_SkippedDOF       = 1 << 1,
    DS_OwnedDOF         = 1 << 2,
    DS_GloballyOwnedDOF = 1 << 3,
    DS_GhostedDOF       = 1 << 4
  };

class TpetraLinearSystem : public LinearSystem
{
public:
  typedef LinSys::GlobalOrdinal GlobalOrdinal;
  typedef LinSys::LocalOrdinal  LocalOrdinal;

  TpetraLinearSystem(
    Realm &realm,
    const unsigned numDof,
    EquationSystem *eqSys,
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
      unsigned numEntities,
      const stk::mesh::Entity* entities,
      const SharedMemView<const double*> & rhs,
      const SharedMemView<const double**> & lhs,
      const SharedMemView<int*> & localIds,
      const SharedMemView<int*> & sortPermutation,
      const char * trace_tag);

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

  /** Reset LHS and RHS for the given set of nodes to 0
   *
   *  @param nodeList A list of STK node entities whose rows are zeroed out
   *  @param beginPos Starting index (usually 0)
   *  @param endPos Terminating index (1 for scalar quantities; nDim for vectors)
   */
  virtual void resetRows(
    const std::vector<stk::mesh::Entity> nodeList,
    const unsigned beginPos,
    const unsigned endPos);

  // Solve
  int solve(stk::mesh::FieldBase * linearSolutionField);
  void loadComplete();
  void writeToFile(const char * filename, bool useOwned=true);
  void printInfo(bool useOwned=true);
  void writeSolutionToFile(const char * filename, bool useOwned=true);
  size_t lookup_myLID(MyLIDMapType& myLIDs, stk::mesh::EntityId entityId, const char* msg=nullptr, stk::mesh::Entity entity = stk::mesh::Entity())
  {
    return myLIDs[entityId];
  }


  int getDofStatus(stk::mesh::Entity node);

private:
  void beginLinearSystemConstruction();

  void checkError(
    const int err_code,
    const char * msg) {}

  void copy_kokkos_unordered_map_to_sorted_vector(const ConnectionSetKK& connectionSetKK,
                                                  ConnectionVec& connectionVec);

  void compute_graph_row_lengths(const ConnectionVec& connectionVec,
                                 LinSys::RowLengths& globallyOwnedRowLengths,
                                 LinSys::RowLengths& locallyOwnedRowLengths);

  void insert_graph_connections(const ConnectionVec& connectionVec,
                                LinSys::Graph& ownedGraph,
                                int ownedOrSharedMask);

  void fill_entity_to_LID_mapping();

  void copy_tpetra_to_stk(
    const Teuchos::RCP<LinSys::Vector> tpetraVector,
    stk::mesh::FieldBase * stkField);

  // This method copies a stk::mesh::field to a tpetra multivector. Each dof/node is written into a different
  // vector in the multivector.
  void copy_stk_to_tpetra(stk::mesh::FieldBase * stkField,
    const Teuchos::RCP<LinSys::MultiVector> tpetraVector);

  int addConnections(const stk::mesh::Entity* entities,const size_t&);
  void expand_unordered_map(unsigned newCapacityNeeded);
  void checkForNaN(bool useOwned);
  bool checkForZeroRow(bool useOwned, bool doThrow, bool doPrint=false);

  ConnectionSetKK connectionSetKK_ ;
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
  LinSys::Matrix::local_matrix_type ownedLocalMatrix_;
  LinSys::Matrix::local_matrix_type globallyOwnedLocalMatrix_;
  host_view_type ownedLocalRhs_;
  host_view_type globallyOwnedLocalRhs_;

  Teuchos::RCP<LinSys::Matrix> globallyOwnedMatrix_;
  Teuchos::RCP<LinSys::Vector> globallyOwnedRhs_;

  Teuchos::RCP<LinSys::Vector> sln_;
  Teuchos::RCP<LinSys::Vector> globalSln_;
  Teuchos::RCP<LinSys::Export> exporter_;
  Teuchos::RCP<LinSys::Import> importer_;

  MyLIDMapType myLIDs_;
  std::vector<LocalOrdinal> entityToLID_;
  LocalOrdinal maxOwnedRowId_; // = num_owned_nodes * numDof_
  LocalOrdinal maxGloballyOwnedRowId_; // = (num_owned_nodes + num_globallyOwned_nodes) * numDof_

  std::vector<int> sortPermutation_;
};

template<typename T1, typename T2>
void copy_kokkos_unordered_map(const Kokkos::UnorderedMap<T1,T2>& src,
                               Kokkos::UnorderedMap<T1,T2>& dest)
{
  if (src.capacity() > dest.capacity()) {
    dest = Kokkos::UnorderedMap<T1,T2>(src.capacity());
  }

  unsigned capacity = src.capacity();
  unsigned fail_count = 0;
  for(unsigned i=0; i<capacity; ++i) {
    if (src.valid_at(i)) {
      auto insert_result = dest.insert(src.key_at(i));
      fail_count += insert_result.failed() ? 1 : 0;
    }
  }
  ThrowRequire(fail_count == 0);
}

int getDofStatus_impl(stk::mesh::Entity node, const Realm& realm);

} // namespace nalu
} // namespace Sierra

#endif
