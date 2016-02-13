/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <EpetraLinearSystem.h>
#include <ContactInfo.h>
#include <ContactManager.h>
#include <HaloInfo.h>
#include <Realm.h>
#include <Simulation.h>
#include <LinearSolver.h>
#include <master_element/MasterElement.h>
#include <NaluEnv.h>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/CPUTime.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_FEVector.h>
#include <Epetra_Vector.h>
//#include <AztecOO.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_VectorOut.h>
#include <EpetraExt_MultiVectorOut.h>

#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <sstream>
#include <limits>

namespace sierra{
namespace nalu{

#define GID_(gid, ndof, idof)  ((ndof)*((gid)-1)+(idof)+1)
#define LID_(lid, ndof, idof)  ((ndof)*((lid))+(idof))

#define DEBUG_EPETRA 0

struct EpetraCompare
{

  EpetraCompare(){}

  bool operator() (const int e0Id, const int e1Id)
  {
    return e0Id < e1Id ;
  }
};

//==========================================================================
// Class Definition
//==========================================================================
// EpetraLinearSystem - hook to Epetra
//==========================================================================
EpetraLinearSystem::EpetraLinearSystem(
  Realm &realm,
  const unsigned numDof,
  const std::string & name,
  LinearSolver * linearSolver)
  : LinearSystem(realm, numDof, name, linearSolver),
    rowMap_(NULL),
    graph_(NULL),
    lhs_(NULL),
    rhs_(NULL),
    sln_(NULL),
    coords_(0)
{
  //throw std::runtime_error("EpetraLinearSystem is deprecated");
}

EpetraLinearSystem::~EpetraLinearSystem()
{
  EpetraLinearSolver *linearSolver = reinterpret_cast<EpetraLinearSolver *>(linearSolver_);
  if(linearSolver->activateML()) linearSolver->destroyML();
  if(linearSolver->activateMueLu()) linearSolver->destroyMueLu();
  delete sln_;
  delete rhs_;
  if (linearSolver->activateMueLu())
    mueLuLHS = Teuchos::null;
  else
    delete lhs_;
  delete graph_;
  delete rowMap_;
  delete[] coords_;
}

void
EpetraLinearSystem::checkError(const int err_code, const char * msg)
{
  if(err_code >= 0) return;
  std::cout << "checkEpetraError(" << err_code << "): " << msg << std::endl;
  ThrowRequire(err_code >= 0);
}

void
EpetraLinearSystem::beginLinearSystemConstruction()
{
  if(inConstruction_) return;
  inConstruction_ = true;
  ThrowRequire(graph_ == NULL);
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const stk::mesh::Selector s_owned = meta_data.locally_owned_part()
      & !stk::mesh::selectUnion(realm_.get_slave_part_vector());
  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_owned );

  int num_owned_nodes = 0;
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ;
        ib != buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    num_owned_nodes += length;
  }
  std::vector<int> owned_rows;
  owned_rows.reserve(num_owned_nodes * numDof_);

  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ;
        ib != buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const stk::mesh::EntityId id = *stk::mesh::field_data(*realm_.naluGlobalId_, b[k]);
      for(unsigned idof=0; idof < numDof_; ++idof) {
        int row = GID_(id , numDof_ , idof);
        owned_rows.push_back(row);
      }
    }
  }

  EpetraLinearSolver *linearSolver = reinterpret_cast<EpetraLinearSolver *>(linearSolver_);
  if (linearSolver->activateMueLu())
      std::sort(owned_rows.begin(), owned_rows.end(), EpetraCompare() );

  const int num_owned_rows = owned_rows.size();
  int num_global_rows = 0;
  stk::all_reduce_sum(bulk_data.parallel(), &num_owned_rows, &num_global_rows, 1);

  const Epetra_MpiComm epetra_comm = bulk_data.parallel();

  rowMap_ = new Epetra_Map(
      num_global_rows,
      num_owned_rows,
      owned_rows.data(),
      1,
      epetra_comm);
  const int approx_cols_per_row = 0;
  bool ignoreNonLocalEntries = false;
  bool buildNonlocalGraph = false;
  graph_ = new Epetra_FECrsGraph(::Copy, *rowMap_,  approx_cols_per_row, ignoreNonLocalEntries, buildNonlocalGraph);


}

void
EpetraLinearSystem::buildNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const stk::mesh::Selector s_owned = meta_data.locally_owned_part()
    & stk::mesh::selectUnion(parts)
    & !stk::mesh::selectUnion(realm_.get_slave_part_vector());
  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_owned );
  int err_code(0);
  std::vector<int> gids(numDof_);
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ;
        ib != buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      stk::mesh::Entity node = b[k];
      const int node_id = *stk::mesh::field_data(*realm_.naluGlobalId_, node);
      for(unsigned d=0; d < numDof_; ++d) {
        gids[d] = GID_(node_id, numDof_ , d);
      }
      // assume full coupling -- all dofs on each node depend on all dofs on each node.
      err_code = graph_->InsertGlobalIndices(gids.size(), gids.data(), gids.size(), gids.data());
      checkError(err_code, "build_edge_to_node_graph/InsertGlobalIndices");
    }
  }
}

void
EpetraLinearSystem::buildEdgeToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const stk::mesh::Selector s_owned = meta_data.locally_owned_part()
        & stk::mesh::selectUnion(parts);
  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::EDGE_RANK, s_owned );
  const int numNodes = 2; // Edges are easy...
  const int numIds = numNodes * numDof_;
  int err_code(0);
  std::vector<int> gids(numIds);
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ;
        ib != buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const unsigned edge_offset = k;
      stk::mesh::Entity const * const edge_nodes = b.begin_nodes(edge_offset);

      // figure out the global dof ids for each dof on each node
      for(int n=0; n < numNodes; ++n) {
        const int node_id = *stk::mesh::field_data(*realm_.naluGlobalId_, edge_nodes[n]);
        for(unsigned d=0; d < numDof_; ++d)
          gids[n*numDof_ + d] = GID_(node_id, numDof_ , d);
      }
      // assume full coupling -- all dofs on each node depend on all dofs on each node.
      err_code = graph_->InsertGlobalIndices(gids.size(), gids.data(), gids.size(), gids.data());
      checkError(err_code, "build_edge_to_node_graph/InsertGlobalIndices");
    }
  }
}

void
EpetraLinearSystem::buildFaceToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const stk::mesh::Selector s_owned = meta_data.locally_owned_part()
        & stk::mesh::selectUnion(parts);
  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( meta_data.side_rank(), s_owned );
  int err_code(0);
  std::vector<int> gids;
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ;
        ib != buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const unsigned face_offset = k;
      stk::mesh::Entity const * face_nodes = b.begin_nodes(face_offset);

      // figure out the global dof ids for each dof on each node
      const int numNodes = b.num_nodes(face_offset);
      const int numIds = numNodes * numDof_;
      gids.resize(numIds);
      for(int n=0; n < numNodes; ++n) {
        const int nodeId = *stk::mesh::field_data(*realm_.naluGlobalId_, face_nodes[n]);
        for(unsigned d=0; d < numDof_; ++d)
          gids[n*numDof_ + d] = GID_(nodeId , numDof_ , d);
      }
      // assume full coupling -- all dofs on each node depend on all dofs on each node.
      err_code = graph_->InsertGlobalIndices(gids.size(), gids.data(), gids.size(), gids.data());
      checkError(err_code, "build_edge_to_node_graph/InsertGlobalIndices");
    }
  }
}

void
EpetraLinearSystem::buildElemToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const stk::mesh::Selector s_owned = meta_data.locally_owned_part()
        & stk::mesh::selectUnion(parts);
  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_owned );
  int err_code(0);
  std::vector<int> gids;
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ;
        ib != buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const unsigned elem_offset = k;
      stk::mesh::Entity const * elem_nodes = b.begin_nodes(elem_offset);

      // figure out the global dof ids for each dof on each node
      const int numNodes = b.num_nodes(elem_offset);
      const int numIds = numNodes * numDof_;
      gids.resize(numIds);
      for(int n=0; n < numNodes; ++n) {
        const int nodeId = *stk::mesh::field_data(*realm_.naluGlobalId_, elem_nodes[n]);
        for(unsigned d=0; d < numDof_; ++d)
          gids[n*numDof_ + d] = GID_(nodeId, numDof_ , d);
      }
      // assume full coupling -- all dofs on each node depend on all dofs on each node.
      err_code = graph_->InsertGlobalIndices(gids.size(), gids.data(), gids.size(), gids.data());
      checkError(err_code, "build_elem_to_node_graph/InsertGlobalIndices");
    }
  }
}

void
EpetraLinearSystem::buildReducedElemToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const stk::mesh::Selector s_owned = meta_data.locally_owned_part()
        & stk::mesh::selectUnion(parts);
  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_owned );
  int err_code(0);
  std::vector<int> gids;
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ;
        ib != buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    // extract master element
    MasterElement *meSCS = realm_.get_surface_master_element(b.topology());
    // extract master element specifics
    const int numScsIp = meSCS->numIntPoints_;
    const int *lrscv = meSCS->adjacentNodes();

    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const unsigned elem_offset = k;
      stk::mesh::Entity const * elem_nodes = b.begin_nodes(elem_offset);

      // figure out the global dof ids for each dof on each node
      const int numIds = 2 * numDof_;
      gids.resize(numIds);
      for(int n=0; n < numScsIp; ++n) {
        int nodeId[2];
        nodeId[0] = *stk::mesh::field_data(*realm_.naluGlobalId_, elem_nodes[lrscv[2*n]]);
        nodeId[1] = *stk::mesh::field_data(*realm_.naluGlobalId_, elem_nodes[lrscv[2*n+1]]);
        for (int j = 0; j < 2; ++j){
          for (unsigned d  = 0; d < numDof_; ++d)
            gids[j*numDof_ + d] = GID_(nodeId[j], numDof_, d);
        }
        // assume full coupling -- all dofs on each node depend on all dofs on each node.
        err_code = graph_->InsertGlobalIndices(gids.size(), gids.data(), gids.size(), gids.data());
        checkError(err_code, "build_reduced_elem_to_node_graph/InsertGlobalIndices");

      }
    }
  }
}


void
EpetraLinearSystem::buildFaceElemToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const stk::mesh::Selector s_owned = meta_data.locally_owned_part()
        & stk::mesh::selectUnion(parts);
  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_owned );
  int err_code(0);
  std::vector<int> gids;
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin() ;
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      const stk::mesh::Entity face = b[k];

      // extract the connected element to this exposed face; should be single in size!
      const stk::mesh::Entity* face_elem_rels = bulk_data.begin_elements(face);//face.relations(stk::topology::ELEMENT_RANK);
      ThrowAssert( bulk_data.num_elements(face) == 1 );

      // get connected element
      stk::mesh::Entity element = face_elem_rels[0];
      const stk::mesh::Entity* elem_nodes = bulk_data.begin_nodes(element);

      // figure out the global dof ids for each dof on each node
      const int numNodes = bulk_data.num_nodes(element);
      const int numIds = numNodes * numDof_;
      gids.resize(numIds);
      for(int n=0; n < numNodes; ++n) {
        const int nodeId = *stk::mesh::field_data(*realm_.naluGlobalId_, elem_nodes[n]);
        for(unsigned d=0; d < numDof_; ++d)
          gids[n*numDof_ + d] = GID_(nodeId, numDof_ , d);
      }
      // assume full coupling -- all dofs on each node depend on all dofs on each node.
      err_code = graph_->InsertGlobalIndices(gids.size(), gids.data(), gids.size(), gids.data());
      checkError(err_code, "build_elem_to_node_graph/InsertGlobalIndices");
    }
  }
}

void
EpetraLinearSystem::buildEdgeHaloNodeGraph(
  const stk::mesh::PartVector &/*parts*/)
{

  beginLinearSystemConstruction();

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  int err_code(0);

  // iterate contactInfoVec_
  std::vector<ContactInfo *>::iterator ii;
  for( ii=realm_.contactManager_->contactInfoVec_.begin();
       ii!=realm_.contactManager_->contactInfoVec_.end(); ++ii ) {

    // iterate halo face nodes
    std::map<uint64_t, HaloInfo *>::iterator iterHalo;
    for (iterHalo  = (*ii)->haloInfoMap_.begin();
         iterHalo != (*ii)->haloInfoMap_.end();
         ++iterHalo) {

      // halo info object of interest
      HaloInfo * infoObject = (*iterHalo).second;

      // extract element mesh object and global id for face node
      stk::mesh::Entity elem = infoObject->owningElement_;
      stk::mesh::Entity node = infoObject->faceNode_;

      // relations

      stk::mesh::Entity const* elem_nodes = bulk_data.begin_nodes(elem);
      const int numNodes = bulk_data.num_nodes(elem);

      const int numIds = (numNodes+1) * numDof_;
      std::vector<int> gids(numIds);

      // figure out the global dof ids for each dof on each node; n = 0
      const int faceNodeId = *stk::mesh::field_data(*realm_.naluGlobalId_, node);
      for(unsigned d=0; d < numDof_; ++d)
        gids[d] = GID_(faceNodeId, numDof_ , d);

      // proceed with element nodes
      for(int n=0; n < numNodes; ++n) {
        const int elemNodeId = *stk::mesh::field_data(*realm_.naluGlobalId_, elem_nodes[n]);
        for(unsigned d=0; d < numDof_; ++d)
          gids[(n+1)*numDof_ + d] = GID_(elemNodeId, numDof_ , d);
      }
      // assume full coupling -- all dofs on each node depend on all dofs on each node.
      err_code = graph_->InsertGlobalIndices(gids.size(), gids.data(), gids.size(), gids.data());
      checkError(err_code, "build_edge_to_node_graph/InsertGlobalIndices");
    }
  }
}

void
EpetraLinearSystem::buildNonConformalNodeGraph(
  const stk::mesh::PartVector &/*parts*/)
{
  throw std::runtime_error("EpetraLinearSystem can not be used with non conformal algorithm; activate tpetra");
}

void
EpetraLinearSystem::buildOversetNodeGraph(
  const stk::mesh::PartVector &/*parts*/)
{
  throw std::runtime_error("EpetraLinearSystem can not be used with overset algorithm; activate tpetra");
}

/// Copy values from Field X to DataVector u.
void copy(
  stk::mesh::BulkData & bulkData,
  GlobalIdFieldType *naluGlobalId,
  stk::mesh::BucketVector const &node_buckets,
  VectorFieldType * X,
  Epetra_MultiVector &U, const int ndofs)
{

  ThrowRequire(ndofs == static_cast<int>(U.NumVectors()) );

  stk::mesh::BucketVector::const_iterator ib     = node_buckets.begin();
  stk::mesh::BucketVector::const_iterator ib_end = node_buckets.end();
  for ( ; ib != ib_end ; ++ib ) {
    const stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length = b.size() ;
    stk::mesh::Bucket::iterator nodes = b.begin();
    for ( stk::mesh::Bucket::size_type n=0; n < length; ++n) {
      stk::mesh::Entity node = nodes[n];
      double * value = stk::mesh::field_data(*X ,node);
      const stk::mesh::EntityId nodeId = *stk::mesh::field_data(*naluGlobalId, node);
      const int nodeIntId = nodeId;
      ThrowRequire(nodeId>0);
      ThrowRequire(NaluEnv::self().parallel_rank() == static_cast<int>(bulkData.parallel_owner_rank(node)));
      for (int ndof=0; ndof<ndofs; ++ndof) {
        U.ReplaceGlobalValue(nodeIntId, ndof, value[ndof]);
      }
    }
  }
}


void
EpetraLinearSystem::finalizeLinearSystem()
{
  ThrowRequire(inConstruction_);
  int err_code(0);
  err_code = graph_->GlobalAssemble();
  checkError(err_code, "finalize_linear_system - FillComplete");
  lhs_ = new Epetra_FECrsMatrix(::Copy, *graph_);
  rhs_ = new Epetra_FEVector(*rowMap_);
  sln_ = new Epetra_Vector(*rowMap_);
  err_code = lhs_->OptimizeStorage();
  checkError(err_code, "finalize_linear_system - OptimizeStorage");
  zeroSystem();
  inConstruction_ = false;

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  EpetraLinearSolver *linearSolver = reinterpret_cast<EpetraLinearSolver *>(linearSolver_);

  if(linearSolver->activateMueLu())
  {
    const int nDim = meta_data.spatial_dimension();

    Teuchos::RCP<Epetra_MultiVector> coords
      = Teuchos::RCP<Epetra_MultiVector >(
      new Epetra_MultiVector(*rowMap_, nDim));

    stk::mesh::Selector s_locally_owned = meta_data.locally_owned_part();

    stk::mesh::BucketVector const & node_buckets =
        bulk_data.get_buckets(stk::topology::NODE_RANK, s_locally_owned);

    VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

    copy(bulk_data, realm_.naluGlobalId_, node_buckets, coordinates, *coords, nDim);

    linearSolver->setMueLuCoordinates(coords);

    mueLuLHS = Teuchos::rcp(lhs_);
    linearSolver->setMueLuMatrix(mueLuLHS);
  }

  if(linearSolver->activateML())
  {
    VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

    // setup for buckets; union parts and ask for locally owned
    stk::mesh::Selector s_locally_owned = meta_data.locally_owned_part();
    stk::mesh::BucketVector const& node_buckets =
      realm_.get_buckets( stk::topology::NODE_RANK, s_locally_owned );

    int num_owned_nodes = 0;
    for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
      ib != node_buckets.end() ; ++ib ) {
      num_owned_nodes += (*ib)->size();
    }

    const int nDim = meta_data.spatial_dimension();

    coords_ = new double[nDim*num_owned_nodes];

    double * xcoords = &coords_[0];
    double * ycoords = nDim > 1 ? &coords_[num_owned_nodes] : NULL;
    double * zcoords = nDim > 2 ? &coords_[2*num_owned_nodes] : NULL;

    stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
    stk::mesh::BucketVector::const_iterator ib_end = node_buckets.end();

    int offset[nDim];
    for (int iDir = 0; iDir < nDim; ++iDir)
      offset[iDir] = iDir*num_owned_nodes;
    // separate coordinates into x, y, and z in same order as is passed to Epetra_Map
    for ( ; ib != ib_end; ++ib )
    {
      const stk::mesh::Bucket & b = **ib;
      const stk::mesh::Bucket::size_type length = b.size();
      stk::mesh::Bucket::iterator nodes = b.begin();
      for (stk::mesh::Bucket::size_type n=0; n < length; ++n)
      {
        stk::mesh::Entity node = nodes[n];
        double * value = stk::mesh::field_data(*coordinates ,node);
        ThrowRequire(*stk::mesh::field_data(*realm_.naluGlobalId_, node) > 0);
        ThrowRequire(NaluEnv::self().parallel_rank() == static_cast<int>(bulk_data.parallel_owner_rank(node)));
        for (int iDir = 0; iDir < nDim; ++iDir)
        {
          coords_[offset[iDir]] = value[iDir];
          ++offset[iDir];
        }

      }
    }

    linearSolver->populateMLCoordinates(xcoords, ycoords, zcoords);

  }
}

void
EpetraLinearSystem::zeroSystem()
{
  int err_code(0);
  err_code = lhs_->PutScalar(0);
  checkError(err_code, "zero - lhs_->PutScalar");

  err_code = rhs_->PutScalar(0);
  checkError(err_code, "zero - rhs_->PutScalar");

  err_code = sln_->PutScalar(0);
  checkError(err_code, "zero - sln_->PutScalar");
}

void
EpetraLinearSystem::sumInto(
  const std::vector<stk::mesh::Entity> & sym_meshobj,
  std::vector<int> &scratchIds,
  std::vector<double> &/*scratchVals*/,
  const std::vector<double> & rhs,
  const std::vector<double> & lhs,
  const char *trace_tag
  )
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  const unsigned p_rank = bulk_data.parallel_rank();
  (void)p_rank;

  int err_code(0);
  const int n_obj = sym_meshobj.size();
  const int numRows = n_obj * numDof_;

  ThrowAssert((size_t)numRows == rhs.size());
  ThrowAssert((size_t)numRows*numRows == lhs.size());

  for(int i=0; i < n_obj; ++i) {
    const int globalOffset = GID_(*stk::mesh::field_data(*realm_.naluGlobalId_, sym_meshobj[i]), numDof_, 0);
    for(unsigned d=0; d < numDof_; ++d) {
      int lid = i*numDof_ + d;
      scratchIds[lid] = globalOffset + d;
    }
  }

  err_code = rhs_->SumIntoGlobalValues(numRows, scratchIds.data(), rhs.data());
  checkError(err_code, "sum_into - rhs->SumIntoGlobalValues");

  err_code = lhs_->SumIntoGlobalValues(numRows, scratchIds.data(), lhs.data(), Epetra_FECrsMatrix::ROW_MAJOR);
  checkError(err_code, "sum_into - lhs->SumIntoGlobalValues");
}

void
EpetraLinearSystem::applyDirichletBCs(
  stk::mesh::FieldBase * solutionField,
  stk::mesh::FieldBase * bcValuesField,
  const stk::mesh::PartVector & parts,
  const unsigned beginPos,
  const unsigned endPos)
{
  double adbc_time = -stk::cpu_time();

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  const unsigned p_size = bulk_data.parallel_size();
  const unsigned p_rank = bulk_data.parallel_rank();
  (void)p_rank;
  (void)p_size;
  int err_code = 0;

  const stk::mesh::Selector selector = stk::mesh::selectUnion(parts) &
    stk::mesh::selectField(*solutionField)
#if 0
    & !stk::mesh::selectUnion(realm_.get_slave_part_vector())
#endif
;

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, selector );

  int nbc=0;
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin();
        ib != buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const unsigned fieldSize = field_bytes_per_entity(*solutionField, b) / sizeof(double);
    ThrowRequire(fieldSize == numDof_);

    bool ghost_bucket = !b.owned() && !b.shared();
    if (ghost_bucket)
      continue;

    const stk::mesh::Bucket::size_type length   = b.size();
    const double * solution = (double*)stk::mesh::field_data(*solutionField, b);
    const double * bcValues = (double*)stk::mesh::field_data(*bcValuesField, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const stk::mesh::Entity entity = b[k];
      const int entity_id = *stk::mesh::field_data(*realm_.naluGlobalId_, entity);

      // FIXME: what about multiple solution fields, e.g., velocity AND pressure
      for(unsigned d=beginPos; d < endPos; ++d) {
        const int rowGID = GID_(entity_id , numDof_ , d);
        const int rowLID = lhs_->Map().LID(rowGID);

        // Adjust the LHS
        int numEntries = lhs_->ColMap().NumMyPoints();
        std::vector<int> colGIDs(numEntries);
        std::vector<double> newVals(numEntries);
        for(int entry=0; entry < numEntries; ++entry)
          {
            const int colGID = lhs_->ColMap().GID(entry);
            colGIDs[entry] =colGID;
            double newv = colGID == rowGID ? 1 : 0;
            if (rowLID < 0) newv = 0;
            newVals[entry] = newv;
          }

        //const int colLID = lhs_->ColMap().LID(colGID);
        err_code = lhs_->ReplaceGlobalValues(rowGID, numEntries, newVals.data(), colGIDs.data());
        checkError(err_code, "LinearSystem::applyDirichletBCs/modify LHS");

        // Replace the RHS residual with (desired - actual)
        {
          ++nbc;
          const double bc_residual = (rowLID < 0 ? 0 : (bcValues[k*fieldSize + d] - solution[k*fieldSize + d]));
          err_code = rhs_->ReplaceGlobalValues(1,  &rowGID, &bc_residual);
          checkError(err_code, "LinearSystem::applyDirichletBCs/modify RHS");
        }
      }
    }
  }

  adbc_time += stk::cpu_time();
}

void
EpetraLinearSystem::loadComplete()
{
  bool reuse_map_and_exporter = false;
  bool callFillComplete = true;
  Epetra_CombineMode mode = Add;
  int err_code=0;
  err_code = lhs_->GlobalAssemble(callFillComplete, mode, reuse_map_and_exporter);
  checkError(err_code, "solve - GlobalAssemble LHS");

  err_code = rhs_->GlobalAssemble(mode, reuse_map_and_exporter);
  checkError(err_code, "solve - GlobalAssemble RHS");

}

int
EpetraLinearSystem::solve(stk::mesh::FieldBase * linearSolutionField)
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  const unsigned p_size = bulk_data.parallel_size();
  const unsigned p_rank = bulk_data.parallel_rank();
  (void)p_rank;
  (void)p_size;

  EpetraLinearSolver *linearSolver = reinterpret_cast<EpetraLinearSolver *>(linearSolver_);
  linearSolver->setSystemObjects(lhs_, rhs_);
  if (linearSolver->getConfig()->getWriteMatrixFiles())
    writeToFile(this->name_.c_str());

  double solve_time = -stk::cpu_time();

  int iters;
  double finalResidNorm;

  // memory diagnostic
  if ( realm_.get_activate_memory_diagnostic() ) {
    NaluEnv::self().naluOutputP0() << "NaluMemory::EpetraLinearSystem::solve() PreSolve: " << name_ << std::endl;
    realm_.provide_memory_summary();
  }

  const int status = linearSolver->solve(
      sln_,
      iters,
      finalResidNorm);
  solve_time += stk::cpu_time();

  if (linearSolver->getConfig()->getWriteMatrixFiles())
  {
    writeSolutionToFile(this->name_.c_str());
    ++writeCounter_;
  }

  copy_epetra_to_stk(sln_, linearSolutionField);
  sync_field(linearSolutionField);

  // computeL2 norm
  std::vector<double> theNorm(1); // for now, size is unity..
  rhs_->Norm2(&theNorm[0]);

  // save off solver info
  linearSolveIterations_ = iters;
  nonLinearResidual_ = realm_.l2Scaling_*theNorm[0];
  linearResidual_ = finalResidNorm;

  if ( realm_.currentNonlinearIteration_ == 1 )
    firstNonLinearResidual_ = nonLinearResidual_;
  scaledNonLinearResidual_ = nonLinearResidual_/std::max(std::numeric_limits<double>::epsilon(), firstNonLinearResidual_);

  if ( provideOutput_ ) {
    const int nameOffset = name_.length()+8;
    NaluEnv::self().naluOutputP0()
      << std::setw(nameOffset) << std::right << name_
      << std::setw(32-nameOffset)  << std::right << iters
      << std::setw(18) << std::right << finalResidNorm
      << std::setw(15) << std::right << nonLinearResidual_
      << std::setw(14) << std::right << scaledNonLinearResidual_ << std::endl;
  }

  return status;
}

void
EpetraLinearSystem::writeToFile(const char * base_filename, bool useOwned)
{
  const int currentCount = writeCounter_;
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  const unsigned p_rank = bulkData.parallel_rank();
  const unsigned p_size = bulkData.parallel_size();

  if (1)
    {
      std::ostringstream osLhs;
      std::ostringstream osRhs;

      osLhs << base_filename << "-" << currentCount << ".epetra.mm." << p_size; // A little hacky but whatever
      int err_code = EpetraExt::RowMatrixToMatrixMarketFile(osLhs.str().c_str(), *lhs_, 0, 0, true);
      checkError(err_code, "RowMatrixToMatrixMarketFile");

      osRhs << base_filename << "-" << currentCount << ".epetra.rhs." << p_size; // A little hacky but whatever
      err_code = EpetraExt::MultiVectorToMatrixMarketFile(osRhs.str().c_str(), *rhs_, 0, 0, true);
      checkError(err_code, "VectorToMatrixMarketFile");
    }

  if (1)
    {
      std::ostringstream osLhs;
      std::ostringstream osRhs;
      std::ostringstream osSln;
      std::ostringstream osGra;

      osLhs << base_filename << "-" << currentCount << ".epetra.mm." << p_size << "." << p_rank; // A little hacky but whatever
      osRhs << base_filename << "-" << currentCount << ".epetra.rhs." << p_size << "." << p_rank; // A little hacky but whatever
      osGra << base_filename << "-" << currentCount << ".epetra.gra." << p_size << "." << p_rank; // A little hacky but whatever
      //osSln << base_filename << "-" << "O-" << currentCount << ".sln." << p_size << "." << p_rank; // A little hacky but whatever
      
    }

}

void
EpetraLinearSystem::writeSolutionToFile(const char * base_filename, bool useOwned)
{
  const int currentCount = writeCounter_;
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  const unsigned p_rank = bulkData.parallel_rank();
  const unsigned p_size = bulkData.parallel_size();

  if (1)
    {
      std::ostringstream osSln;

      osSln << base_filename << "-" << currentCount << ".epetra.sln." << p_size; // A little hacky but whatever
      int err_code = EpetraExt::VectorToMatrixMarketFile(osSln.str().c_str(), *sln_, 0, 0, true);
      checkError(err_code, "VectorToMatrixMarketFile");
    }

  if (1)
    {
      std::ostringstream osSln;

      osSln << base_filename << "-" << currentCount << ".epetra.sln." << p_size << "." << p_rank; // A little hacky but whatever

    }

}

void
EpetraLinearSystem::copy_epetra_to_stk(
  const Epetra_Vector * epetraField,
  stk::mesh::FieldBase * stkField)
{
  ThrowAssert(epetraField);
  ThrowAssert(stkField);
  const Epetra_Vector & epetraVector = *epetraField;

  const stk::mesh::Selector selector = stk::mesh::selectField(*stkField)
    & !stk::mesh::selectUnion(realm_.get_slave_part_vector());

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets(stk::topology::NODE_RANK, selector);

  for (size_t ib=0; ib < buckets.size(); ++ib) {
    stk::mesh::Bucket & b = *buckets[ib];

    const unsigned fieldSize = field_bytes_per_entity(*stkField, b) / sizeof(double);
    ThrowRequire(fieldSize == numDof_);

    if (b.owned())   {
      const stk::mesh::Bucket::size_type length = b.size();
      double * stkFieldPtr = (double*)stk::mesh::field_data(*stkField, b);
      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        const stk::mesh::Entity entity = b[k];
        const int entity_id = *stk::mesh::field_data(*realm_.naluGlobalId_, entity);
        for(unsigned d=0; d < fieldSize; ++d) {
          const int epetraGID = GID_(entity_id , numDof_ , d);
          const int epetraLID = epetraVector.Map().LID(epetraGID);
          ThrowRequire(epetraLID >= 0);
          const int stkIndex = k*numDof_ + d;
          stkFieldPtr[stkIndex] = epetraVector[epetraLID];
        }
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
