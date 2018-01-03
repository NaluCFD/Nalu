/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <TpetraLinearSystem.h>
#include <NonConformalInfo.h>
#include <NonConformalManager.h>
#include <FieldTypeDef.h>
#include <DgInfo.h>
#include <Realm.h>
#include <PeriodicManager.h>
#include <Simulation.h>
#include <LinearSolver.h>
#include <master_element/MasterElement.h>
#include <EquationSystem.h>
#include <NaluEnv.h>

#include <KokkosInterface.h>

// overset
#include <overset/OversetManager.h>
#include <overset/OversetInfo.h>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/WallTime.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

// For Tpetra support
#include <Kokkos_Serial.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Details_shortSort.hpp>

#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Tpetra_MatrixIO.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <set>
#include <limits>
#include <type_traits>

#include <sstream>
#define KK_MAP
namespace sierra{
namespace nalu{

#define GID_(gid, ndof, idof)  ((ndof)*((gid)-1)+(idof)+1)
#define LID_(lid, ndof, idof)  ((ndof)*((lid))+(idof))

#define GLOBAL_ENTITY_ID(gid, ndof) ((gid-1)/ndof + 1)
#define GLOBAL_ENTITY_ID_IDOF(gid, ndof) ((gid-1) % ndof)

///====================================================================================================================================
///======== T P E T R A ===============================================================================================================
///====================================================================================================================================

//==========================================================================
// Class Definition
//==========================================================================
// TpetraLinearSystem - hook to Tpetra
//==========================================================================
TpetraLinearSystem::TpetraLinearSystem(
  Realm &realm,
  const unsigned numDof,
  EquationSystem *eqSys,
  LinearSolver * linearSolver)
  : LinearSystem(realm, numDof, eqSys, linearSolver)
{
  Teuchos::ParameterList junk;
  node_ = Teuchos::rcp(new LinSys::Node(junk));
}

TpetraLinearSystem::~TpetraLinearSystem()
{
  // dereference linear solver in safe manner
  TpetraLinearSolver *linearSolver = reinterpret_cast<TpetraLinearSolver *>(linearSolver_);
  linearSolver->destroyLinearSolver();
}

struct CompareEntityById
{
  const stk::mesh::BulkData &m_mesh;
  const GlobalIdFieldType *m_naluGlobalId;

  CompareEntityById(
    const stk::mesh::BulkData &mesh, const GlobalIdFieldType *naluGlobalId)
    : m_mesh(mesh),
      m_naluGlobalId(naluGlobalId) {}

  bool operator() (const stk::mesh::Entity& e0, const stk::mesh::Entity& e1)
  {
    const stk::mesh::EntityId e0Id = *stk::mesh::field_data(*m_naluGlobalId, e0);
    const stk::mesh::EntityId e1Id = *stk::mesh::field_data(*m_naluGlobalId, e1);
    return e0Id < e1Id ;
  }
  bool operator() (const Connection& c0, const Connection& c1)
  {
    const stk::mesh::EntityId c0firstId = *stk::mesh::field_data(*m_naluGlobalId, c0.first);
    const stk::mesh::EntityId c1firstId = *stk::mesh::field_data(*m_naluGlobalId, c1.first);
    if (c0firstId != c1firstId) {
      return c0firstId < c1firstId;
    }
    const stk::mesh::EntityId c0secondId = *stk::mesh::field_data(*m_naluGlobalId, c0.second);
    const stk::mesh::EntityId c1secondId = *stk::mesh::field_data(*m_naluGlobalId, c1.second);
    return c0secondId < c1secondId;
  }
};

// determines whether the node is to be put into which map/graph/matrix
// FIXME - note that the DOFStatus enum can be Or'd together if need be to
//   distinguish ever more complicated situations, for example, a DOF that
//   is both owned and ghosted: OwnedDOF | GhostedDOF
int TpetraLinearSystem::getDofStatus(stk::mesh::Entity node)
{
    return getDofStatus_impl(node, realm_);
}

void
TpetraLinearSystem::beginLinearSystemConstruction()
{
  if(inConstruction_) return;
  inConstruction_ = true;
  ThrowRequire(ownedGraph_.is_null());
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  stk::mesh::MetaData & metaData = realm_.meta_data();
  const unsigned p_rank = bulkData.parallel_rank();
  (void)p_rank;

  // create a localID for all active nodes in the mesh...
  const stk::mesh::Selector s_universal = metaData.universal_part()
      & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& buckets =
      realm_.get_buckets( stk::topology::NODE_RANK, s_universal );

  // we allow for ghosted nodes when nonconformal is active. When periodic is active, we may
  // also have ghosted nodes due to the periodicGhosting. However, we want to exclude these
  // nodes

  LocalOrdinal numGhostNodes = 0;
  LocalOrdinal numOwnedNodes = 0;
  LocalOrdinal numNodes = 0;
  LocalOrdinal numGloballyOwnedNotLocallyOwned = 0; // these are nodes on other procs
  // First, get the number of owned and globallyOwned (or num_globallyOwned_nodes = num_nodes - num_owned_nodes)
  //KOKKOS: BucketLoop parallel "reduce" is accumulating 4 sums
  kokkos_parallel_for("Nalu::TpetraLinearSystem::beginLinearSystemConstructionA", buckets.size(), [&] (const int& ib) {
    stk::mesh::Bucket & b = *buckets[ib];
    const stk::mesh::Bucket::size_type length = b.size();
    //KOKKOS: intra BucketLoop parallel reduce
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get node
      stk::mesh::Entity node = b[k];
      int status = getDofStatus(node);

      if (status & DS_SkippedDOF)
        continue;

      if (status & DS_OwnedDOF) {
        numNodes++;
        numOwnedNodes++;
      }

      if (status & DS_GloballyOwnedDOF) {
        numNodes++;
        numGloballyOwnedNotLocallyOwned++;
      }

      if (status & DS_GhostedDOF) {
        numGhostNodes++;
      }
    }
  });

  maxOwnedRowId_ = numOwnedNodes * numDof_;
  maxGloballyOwnedRowId_ = numNodes * numDof_;

  // Next, grab all the global ids, owned first, then globallyOwned.
  totalGids_.clear();
  totalGids_.reserve(numNodes * numDof_);

  // Also, we'll build up our own local id map. Note: first we number
  // the owned nodes then we number the globallyOwned nodes.
  LocalOrdinal localId = 0;

  // make separate arrays that hold the owned and globallyOwned gids
  std::vector<stk::mesh::Entity> owned_nodes, globally_owned_nodes;
  std::vector<GlobalOrdinal> ownedGids, globallyOwnedGids;

  // owned first:
  for(const stk::mesh::Bucket* bptr : buckets) {
    const stk::mesh::Bucket & b = *bptr;
    for ( stk::mesh::Entity entity : b ) {
      int status = getDofStatus(entity);
      if (!(status & DS_SkippedDOF) && (status & DS_OwnedDOF))
        owned_nodes.push_back(entity);
    }
  }

  std::sort(owned_nodes.begin(), owned_nodes.end(), CompareEntityById(bulkData, realm_.naluGlobalId_) );

  myLIDs_.clear();
  //KOKKOS: Loop noparallel push_back totalGids_ (std::vector)
  for(stk::mesh::Entity entity : owned_nodes) {
    const stk::mesh::EntityId entityId = *stk::mesh::field_data(*realm_.naluGlobalId_, entity);
    // entityId can be duplicated in periodic or nonconformal
    MyLIDMapType::iterator found = myLIDs_.find(entityId);
    if (found == myLIDs_.end()) {
      myLIDs_[entityId] = localId++;
      for(unsigned idof=0; idof < numDof_; ++ idof) {
        const GlobalOrdinal gid = GID_(entityId, numDof_, idof);
        totalGids_.push_back(gid);
        ownedGids.push_back(gid);
      }
    }
  }
  ThrowRequire(localId == numOwnedNodes);
  
  // now globallyOwned:
  for(const stk::mesh::Bucket* bptr : buckets) {
    const stk::mesh::Bucket & b = *bptr;
    for ( stk::mesh::Entity node : b) {
      int status = getDofStatus(node);
      if (!(status & DS_SkippedDOF) && (status & DS_GloballyOwnedDOF))
        globally_owned_nodes.push_back(node);
    }
  }
  std::sort(globally_owned_nodes.begin(), globally_owned_nodes.end(), CompareEntityById(bulkData, realm_.naluGlobalId_) );

  for (unsigned inode=0; inode < globally_owned_nodes.size(); ++inode) {
    const stk::mesh::Entity entity = globally_owned_nodes[inode];
    const stk::mesh::EntityId naluId = *stk::mesh::field_data(*realm_.naluGlobalId_, entity);
    MyLIDMapType::iterator found = myLIDs_.find(naluId);
    if (found == myLIDs_.end()) {
      myLIDs_[naluId] = localId++;
      for(unsigned idof=0; idof < numDof_; ++ idof) {
        const GlobalOrdinal gid = GID_(naluId, numDof_, idof);
        totalGids_.push_back(gid);
        globallyOwnedGids.push_back(gid);
      }
    }
  }
  
  const int numOwnedRows = numOwnedNodes * numDof_;
  (void)numOwnedRows;

  std::sort(ownedGids.begin(), ownedGids.end());
  std::sort(globallyOwnedGids.begin(), globallyOwnedGids.end());

  const Teuchos::RCP<LinSys::Comm> tpetraComm = Teuchos::rcp(new LinSys::Comm(bulkData.parallel()));
  ownedRowsMap_ = Teuchos::rcp(new LinSys::Map(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), ownedGids, 1, tpetraComm, node_));
  globallyOwnedRowsMap_ = Teuchos::rcp(new LinSys::Map(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), globallyOwnedGids, 1, tpetraComm, node_));

  exporter_ = Teuchos::rcp(new LinSys::Export(globallyOwnedRowsMap_, ownedRowsMap_));
  importer_ = Teuchos::rcp(new LinSys::Import(ownedRowsMap_, globallyOwnedRowsMap_));

  ownedPlusGloballyOwnedRowsMap_ = Teuchos::rcp(new LinSys::Map(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), totalGids_, 1, tpetraComm, node_));

  // Now, we're ready to have Algs call the build*Graph() methods and build up the connection list (row,col).
  size_t numGraphEntriesGuess = 15*numOwnedNodes;
  if (eqSys_->num_graph_entries_ > 0) {
    numGraphEntriesGuess = 1.1*eqSys_->num_graph_entries_;
  }

  connectionSetKK_ = ConnectionSetKK(numGraphEntriesGuess);
}

int TpetraLinearSystem::addConnections(const stk::mesh::Entity* entities, const size_t& num_entities)
{
  int fail_count = 0;
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const unsigned p_rank = bulkData.parallel_rank();
  (void)p_rank;

  for(size_t a=0; a < num_entities; ++a) {
    const stk::mesh::Entity entity_a = entities[a];
    const stk::mesh::EntityId id_a = *stk::mesh::field_data(*realm_.naluGlobalId_, entity_a);
    auto insert_result = connectionSetKK_.insert( Connection(entity_a, entity_a) );
    fail_count += insert_result.failed() ? 1 : 0;

    for(size_t b=a+1; b < num_entities; ++b) {
      const stk::mesh::Entity entity_b = entities[b];
      const stk::mesh::EntityId id_b = *stk::mesh::field_data(*realm_.naluGlobalId_, entity_b);
      const bool a_then_b = id_a < id_b;
      const stk::mesh::Entity entity_min = a_then_b ? entity_a : entity_b;
      const stk::mesh::Entity entity_max = a_then_b ? entity_b : entity_a;
      insert_result = connectionSetKK_.insert( Connection(entity_min, entity_max) );
      fail_count += insert_result.failed() ? 1 : 0;
    }
  }
  return fail_count;
}

void TpetraLinearSystem::expand_unordered_map(unsigned newCapacityNeeded)
{
  ConnectionSetKK tmp(connectionSetKK_.capacity()+newCapacityNeeded);
  copy_kokkos_unordered_map(connectionSetKK_, tmp);
  connectionSetKK_ = tmp;
}

void
TpetraLinearSystem::buildNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & metaData = realm_.meta_data();
//if (realm_.bulk_data().parallel_rank()==0) std::cerr<<"buildNodeGraph"<<std::endl;

  const stk::mesh::Selector s_owned = metaData.locally_owned_part()
    & stk::mesh::selectUnion(parts) 
    & !(stk::mesh::selectUnion(realm_.get_slave_part_vector()))
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_owned );
  int insert_failure_count = 0;
  kokkos_parallel_reduce(buckets.size(), [&] (const int& ib, int& fail_count) {
    const stk::mesh::Bucket & b = *buckets[ib];
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      stk::mesh::Entity node = b[k];
      fail_count += addConnections(&node, 1);
    }
  }, insert_failure_count, "Nalu::TpetraLinearSystem::buildNodeGraph");
  if(insert_failure_count > 0)
  {
//std::cerr<<"buildNodeGraph expanding unordered map, proc "<<realm_.bulkData_->parallel_rank()<<std::endl;
    expand_unordered_map(insert_failure_count);
    buildNodeGraph(parts);
  }
}

void
TpetraLinearSystem::buildEdgeToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & metaData = realm_.meta_data();
//if (realm_.bulk_data().parallel_rank()==0) std::cerr<<"buildEdgeToNodeGraph"<<std::endl;

  const stk::mesh::Selector s_owned = metaData.locally_owned_part()
    & stk::mesh::selectUnion(parts) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::EDGE_RANK, s_owned );
  const size_t numNodes = 2; // Edges are easy...
  std::vector<stk::mesh::Entity> entities(numNodes);
  int insert_failure_count = 0;
  kokkos_parallel_reduce(buckets.size(), [&] (const int& ib, int& fail_count) {
    const stk::mesh::Bucket & b = *buckets[ib];
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      stk::mesh::Entity const * edge_nodes = b.begin_nodes(k);

      // figure out the global dof ids for each dof on each node
      //KOKKOS: nested Loop parallel
      for(size_t n=0; n < numNodes; ++n) {
        entities[n] = edge_nodes[n];
      }
      fail_count += addConnections(entities.data(), entities.size());
    }
  }, insert_failure_count, "Nalu::TpetraLinearSystem::buildEdgeToNodeGraph");
  if(insert_failure_count > 0)
  {
//std::cerr<<"buildEdgeToNodeGraph expanding unordered map, proc "<<realm_.bulkData_->parallel_rank()<<std::endl;
    expand_unordered_map(insert_failure_count);
    buildEdgeToNodeGraph(parts);
  }
}

void
TpetraLinearSystem::buildFaceToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & metaData = realm_.meta_data();
//if (realm_.bulk_data().parallel_rank()==0) std::cerr<<"buildFaceToNodeGraph"<<std::endl;

  const stk::mesh::Selector s_owned = metaData.locally_owned_part()
    & stk::mesh::selectUnion(parts) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( metaData.side_rank(), s_owned );
  int insert_failure_count = 0;
  kokkos_parallel_reduce(buckets.size(), [&] (const int& ib, int & fail_count) {
    const stk::mesh::Bucket & b = *buckets[ib];
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      stk::mesh::Entity const * face_nodes = b.begin_nodes(k);

      const size_t numNodes = b.num_nodes(k);
      fail_count += addConnections(face_nodes, numNodes);
    }
  }, insert_failure_count, "Nalu::TpetraLinearSystem::buildFaceToNodeGraph");
  if(insert_failure_count > 0)
  {
//std::cerr<<"buildFaceToNodeGraph expanding unordered map, proc "<<realm_.bulkData_->parallel_rank()<<std::endl;
    expand_unordered_map(insert_failure_count);
    buildFaceToNodeGraph(parts);
  }
}

void
TpetraLinearSystem::buildElemToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & metaData = realm_.meta_data();
//if (realm_.bulk_data().parallel_rank()==0) std::cerr<<"buildElemToNodeGraph"<<std::endl;

  const stk::mesh::Selector s_owned = metaData.locally_owned_part()
    & stk::mesh::selectUnion(parts) 
    & !(realm_.get_inactive_selector());

  connectionSetKK_.reset_failed_insert_flag();
  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_owned );
  int insert_failure_count = 0;
  for(const stk::mesh::Bucket* bptr : buckets) {
    const stk::mesh::Bucket & b = *bptr;
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      stk::mesh::Entity const * elem_nodes = b.begin_nodes(k);
      const size_t numNodes = b.num_nodes(k);
      insert_failure_count += addConnections(elem_nodes,numNodes);
    }
  }
  if(insert_failure_count > 0)
  {
//std::cerr<<"buildElemToNodeGraph expanding unordered map, proc "<<realm_.bulkData_->parallel_rank()<<std::endl;
    expand_unordered_map(insert_failure_count);
    buildElemToNodeGraph(parts);
  }
  //std::cout<<"KK HashTableSize: "<<connectionSetKK_.size() << " failed_inserts:" << connectionSetKK_.failed_insert() << " load_factor:" << 1.0*connectionSetKK_.size()/connectionSetKK_.capacity() << std::endl;
}

void
TpetraLinearSystem::buildReducedElemToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & metaData = realm_.meta_data();
//if (realm_.bulk_data().parallel_rank()==0) std::cerr<<"buildReducedElemToNodeGraph"<<std::endl;

  const stk::mesh::Selector s_owned = metaData.locally_owned_part()
    & stk::mesh::selectUnion(parts) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_owned );
  std::vector<stk::mesh::Entity> entities;
  int insert_failure_count = 0;
  kokkos_parallel_reduce(buckets.size(), [&] (const int& ib, int& fail_count) {
    const stk::mesh::Bucket & b = *buckets[ib];

    // extract master element
    MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());
    // extract master element specifics
    const int numScsIp = meSCS->numIntPoints_;
    const int *lrscv = meSCS->adjacentNodes();

    const stk::mesh::Bucket::size_type length   = b.size();
    //KOKKOS: intra BucketLoop noparallel addConnections insert (std::set)
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      stk::mesh::Entity const * elem_nodes = b.begin_nodes(k);

      // figure out the global dof ids for each dof on each node
      const size_t numNodes = 2;
      entities.resize(numNodes);
      //KOKKOS: nested Loop noparallel addConnections insert (std::set)
      for (int j = 0; j < numScsIp; ++j){
        //KOKKOS: nested Loop parallel
        for(size_t n=0; n < numNodes; ++n) {
          entities[n] = elem_nodes[lrscv[2*j+n]];
        }
        fail_count += addConnections(entities.data(), entities.size());
      }
    }
  }, insert_failure_count, "Nalu::TpetraLinearSystem::buildReducedElemToNodeGraph");
  if(insert_failure_count > 0)
  {
//std::cerr<<"buildReducedElemToNodeGraph expanding unordered map, proc "<<realm_.bulkData_->parallel_rank()<<std::endl;
    expand_unordered_map(insert_failure_count);
    buildReducedElemToNodeGraph(parts);
  }
}

void
TpetraLinearSystem::buildFaceElemToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  stk::mesh::MetaData & metaData = realm_.meta_data();
//if (realm_.bulk_data().parallel_rank()==0) std::cerr<<"buildFaceElemToNodeGraph"<<std::endl;

  const stk::mesh::Selector s_owned = metaData.locally_owned_part()
    & stk::mesh::selectUnion(parts) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( metaData.side_rank(), s_owned );
  std::vector<stk::mesh::Entity> entities;
  int insert_failure_count = 0;
  kokkos_parallel_reduce(face_buckets.size(), [&] (const int& ib, int& fail_count) {
    const stk::mesh::Bucket & b = *face_buckets[ib];
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const stk::mesh::Entity face = b[k];

      // extract the connected element to this exposed face; should be single in size!
      const stk::mesh::Entity* face_elem_rels = bulkData.begin_elements(face);
      ThrowAssert( bulkData.num_elements(face) == 1 );

      // get connected element and nodal relations
      stk::mesh::Entity element = face_elem_rels[0];
      const stk::mesh::Entity* elem_nodes = bulkData.begin_nodes(element);

      // figure out the global dof ids for each dof on each node
      const size_t numNodes = bulkData.num_nodes(element);
      entities.resize(numNodes);
      //KOKKOS: nested Loop parallel
      for(size_t n=0; n < numNodes; ++n) {
        entities[n] = elem_nodes[n];
      }
      fail_count += addConnections(entities.data(), entities.size());
    }
  }, insert_failure_count, "Nalu::TpetraLinearSystem::buildFaceElemToNodeGraph");
  if(insert_failure_count > 0)
  {
//std::cerr<<"buildFaceElemToNodeGraph expanding unordered map, proc "<<realm_.bulkData_->parallel_rank()<<std::endl;
    expand_unordered_map(insert_failure_count);
    buildFaceElemToNodeGraph(parts);
  }
}

void
TpetraLinearSystem::buildNonConformalNodeGraph(const stk::mesh::PartVector &parts)
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  beginLinearSystemConstruction();
//if (realm_.bulk_data().parallel_rank()==0) std::cerr<<"buildNonConformalNodeGraph"<<std::endl;

  std::vector<stk::mesh::Entity> entities;

  int insert_failure_count = 0;
  // iterate nonConformalManager's dgInfoVec
  std::vector<NonConformalInfo *>::iterator ii;
  for( ii=realm_.nonConformalManager_->nonConformalInfoVec_.begin();
       ii!=realm_.nonConformalManager_->nonConformalInfoVec_.end(); ++ii ) {

    // extract vector of DgInfo
    std::vector<std::vector<DgInfo *> > &dgInfoVec = (*ii)->dgInfoVec_;
    
    std::vector<std::vector<DgInfo*> >::iterator idg;
    for( idg=dgInfoVec.begin(); idg!=dgInfoVec.end(); ++idg ) {

      std::vector<DgInfo *> &faceDgInfoVec = (*idg);

      // now loop over all the DgInfo objects on this particular exposed face
      for ( size_t k = 0; k < faceDgInfoVec.size(); ++k ) {

        DgInfo *dgInfo = faceDgInfoVec[k];

        // extract current/opposing element
        stk::mesh::Entity currentElement = dgInfo->currentElement_;
        stk::mesh::Entity opposingElement = dgInfo->opposingElement_;
        
        // node relations; current and opposing
        stk::mesh::Entity const* current_elem_node_rels = bulkData.begin_nodes(currentElement);
        const int current_num_elem_nodes = bulkData.num_nodes(currentElement);
        stk::mesh::Entity const* opposing_elem_node_rels = bulkData.begin_nodes(opposingElement);
        const int opposing_num_elem_nodes = bulkData.num_nodes(opposingElement);
        
        // resize based on both current and opposing face node size
        entities.resize(current_num_elem_nodes+opposing_num_elem_nodes);
        
        // fill in connected nodes; current
        //KOKKOS: nested Loop parallel
        for ( int ni = 0; ni < current_num_elem_nodes; ++ni ) {
          entities[ni] = current_elem_node_rels[ni];
        }
        
        // fill in connected nodes; opposing
        //KOKKOS: nested Loop parallel
        for ( int ni = 0; ni < opposing_num_elem_nodes; ++ni ) {
          entities[current_num_elem_nodes+ni] = opposing_elem_node_rels[ni];
        }
        
        // okay, now add the connections; will be symmetric 
        // columns of current node row (opposing nodes) will add columns to opposing nodes row
        insert_failure_count += addConnections(entities.data(), entities.size());
      }
    }
  }
  if (insert_failure_count > 0) {
//std::cerr<<"buildNonConformalNodeGraph expanding unordered map, proc "<<realm_.bulkData_->parallel_rank()<<std::endl;
    expand_unordered_map(insert_failure_count);
    buildNonConformalNodeGraph(parts);
  }
}

void
TpetraLinearSystem::buildOversetNodeGraph(const stk::mesh::PartVector &parts)
{
  // extract the rank
  const int theRank = NaluEnv::self().parallel_rank();

  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  beginLinearSystemConstruction();
//if (realm_.bulk_data().parallel_rank()==0) std::cerr<<"buildOversetNodeGraph"<<std::endl;

  std::vector<stk::mesh::Entity> entities;

  int insert_failure_count = 0;
  // iterate oversetInfoVec_
  std::vector<OversetInfo *>::iterator ii;
  for( ii=realm_.oversetManager_->oversetInfoVec_.begin();
       ii!=realm_.oversetManager_->oversetInfoVec_.end(); ++ii ) {

    // extract element mesh object and orphan node
    stk::mesh::Entity owningElement = (*ii)->owningElement_;
    stk::mesh::Entity orphanNode = (*ii)->orphanNode_;

    // extract the owning rank for this node
    const int nodeRank = bulkData.parallel_owner_rank(orphanNode);

    // check to see if this node is locally owned by this rank; we only want to process locally owned nodes, not shared
    if ( theRank != nodeRank )
      continue;

    // relations
    stk::mesh::Entity const* elem_nodes = bulkData.begin_nodes(owningElement);
    const size_t numNodes = bulkData.num_nodes(owningElement);
    const size_t numEntities = numNodes+1;
    entities.resize(numEntities);
    
    entities[0] = orphanNode;
    for(size_t n=0; n < numNodes; ++n) {
      entities[n+1] = elem_nodes[n];
    }
    insert_failure_count += addConnections(entities.data(), entities.size());
  }
  if (insert_failure_count > 0) {
//std::cerr<<"buildOversetNodeGraph expanding unordered map, proc "<<realm_.bulkData_->parallel_rank()<<std::endl;
    expand_unordered_map(insert_failure_count);
    buildOversetNodeGraph(parts);
  }
}

void
TpetraLinearSystem::copy_stk_to_tpetra(
  stk::mesh::FieldBase * stkField,
  const Teuchos::RCP<LinSys::MultiVector> tpetraField)
{
  ThrowAssert(!tpetraField.is_null());
  ThrowAssert(stkField);
  const int numVectors = tpetraField->getNumVectors();

  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  stk::mesh::MetaData & metaData = realm_.meta_data();

  const stk::mesh::Selector selector = stk::mesh::selectField(*stkField) 
    & metaData.locally_owned_part() 
    & !(stk::mesh::selectUnion(realm_.get_slave_part_vector())) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& buckets = bulkData.get_buckets(stk::topology::NODE_RANK, selector);

  for(const stk::mesh::Bucket* bptr : buckets) {
    const stk::mesh::Bucket & b = *bptr;

    const int fieldSize = field_bytes_per_entity(*stkField, b) / (sizeof(double));

    ThrowRequire(numVectors == fieldSize);

    const stk::mesh::Bucket::size_type length = b.size();

    const double * stkFieldPtr = (double*)stk::mesh::field_data(*stkField, b);

    for (stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k )
    {
      const stk::mesh::Entity node = b[k];

      int status = getDofStatus(node);
      if ((status & DS_SkippedDOF) || (status & DS_GloballyOwnedDOF))
        continue;

      const stk::mesh::EntityId nodeId = *stk::mesh::field_data(*realm_.naluGlobalId_, node);
      for(int d=0; d < fieldSize; ++d)
      {
        const size_t stkIndex = k*fieldSize + d;
        tpetraField->replaceGlobalValue(nodeId, d, stkFieldPtr[stkIndex]);
      }
    }
  }
}

void
TpetraLinearSystem::copy_kokkos_unordered_map_to_sorted_vector(const ConnectionSetKK& connectionSetKK,
                                                               ConnectionVec& connectionVec)
{
  connectionVec.resize(connectionSetKK.size());
  Kokkos::parallel_scan("Nalu::TpetraLinearSystem::copy_kokkos_unordered_map_to_sorted_vector",
      connectionSetKK.capacity(), [&] (const int& i, int & ikk, const bool final) {
    if(connectionSetKK.valid_at(i))
    {
      if(final)
        connectionVec[ikk] = connectionSetKK.key_at(i);

      ikk++;
    }
  });
  size_t actualNumGraphEntries = connectionSetKK.size();
  //std::cerr<<"finalize, actual graph entries: "<<connectionSetKK_.size()<<std::endl;
  if (actualNumGraphEntries > eqSys_->num_graph_entries_) {
    eqSys_->num_graph_entries_ = actualNumGraphEntries;
  }

  std::sort(connectionVec.begin(), connectionVec.end());
}

void
TpetraLinearSystem::compute_graph_row_lengths(const ConnectionVec& connectionVec,
                                              LinSys::RowLengths& globallyOwnedRowLengths,
                                              LinSys::RowLengths& locallyOwnedRowLengths)
{
  Kokkos::View<size_t*,DeviceSpace> deviceGloballyOwnedRowLengths = globallyOwnedRowLengths.view<DeviceSpace>();
  Kokkos::View<size_t*,DeviceSpace> deviceLocallyOwnedRowLengths = locallyOwnedRowLengths.view<DeviceSpace>();

  size_t numConnections = connectionVec.size();
  kokkos_parallel_for("Nalu::TpetraLinearSystem::compute_graph_row_lengths", numConnections, [&] (const size_t& i)
  {
    const stk::mesh::Entity entity_a = connectionVec[i].first;
    const stk::mesh::Entity entity_b = connectionVec[i].second;

    const stk::mesh::EntityId entityId_a = *stk::mesh::field_data(*realm_.naluGlobalId_, entity_a);
    const stk::mesh::EntityId entityId_b = *stk::mesh::field_data(*realm_.naluGlobalId_, entity_b);

    const int entity_a_status = getDofStatus(entity_a);
    if (entity_a_status & DS_GloballyOwnedDOF) { // !Locally owned
      LocalOrdinal lid_a = globallyOwnedRowsMap_->getLocalElement(GID_(entityId_a, numDof_ , 0));
      const bool add_b_to_a_connection = (entity_a != entity_b) && (getDofStatus(entity_b) & DS_GloballyOwnedDOF);
      LocalOrdinal lid_b = add_b_to_a_connection ? globallyOwnedRowsMap_->getLocalElement(GID_(entityId_b, numDof_, 0)) : 0;

      for (size_t d=0; d < numDof_; ++d) {
        deviceGloballyOwnedRowLengths(lid_a+d) += numDof_;
  
        if (add_b_to_a_connection) {
          deviceGloballyOwnedRowLengths(lid_b+d) += numDof_;
        }
      }
    }

    if (entity_a_status & DS_OwnedDOF) { // Locally owned
      LocalOrdinal lid_a = ownedRowsMap_->getLocalElement(GID_(entityId_a, numDof_ , 0));
      const bool add_b_to_a_connection = (entity_a != entity_b) && (getDofStatus(entity_b) & DS_OwnedDOF);
      LocalOrdinal lid_b = add_b_to_a_connection ? ownedRowsMap_->getLocalElement(GID_(entityId_b, numDof_, 0)) : 0;

      for (size_t d=0; d < numDof_; ++d) {
        deviceLocallyOwnedRowLengths(lid_a+d) += numDof_;
  
        if (add_b_to_a_connection) {
          deviceLocallyOwnedRowLengths(lid_b+d) += numDof_;
        }
      }
    }
  });
}

void
TpetraLinearSystem::insert_graph_connections(const ConnectionVec& connectionVec,
                                             LinSys::Graph& crsGraph,
                                             int ownedOrSharedMask)
{
  std::vector<GlobalOrdinal> globalDofs_a(numDof_);
  const unsigned max = 64;
  std::vector<stk::mesh::Entity> entities_b(max);
  std::vector<stk::mesh::EntityId> entityIds_b(max);
  std::vector<GlobalOrdinal> globalDofs_b(max*numDof_);
  const size_t numConnections = connectionVec.size();

  unsigned numColEntities = 0;
  // Insert all the local connection data
  //KOKKOS: Loop noparallel Graph insertGlobalIndices
  for(size_t i=0; i<numConnections; ) {
    globalDofs_b.resize(max*numDof_);
    const stk::mesh::Entity entity_a = connectionVec[i].first;
    const stk::mesh::EntityId entityId_a = *stk::mesh::field_data(*realm_.naluGlobalId_, entity_a);
    for (size_t d=0; d < numDof_; ++d) {
      globalDofs_a[d] = GID_(entityId_a, numDof_, d);
    }

    while(numColEntities<max && i+numColEntities<numConnections && connectionVec[i+numColEntities].first == entity_a) {
        const stk::mesh::Entity entity_b = connectionVec[i+numColEntities].second;
        entities_b[numColEntities] = entity_b;
        const stk::mesh::EntityId entityId_b = *stk::mesh::field_data(*realm_.naluGlobalId_, entity_b);
        entityIds_b[numColEntities] = entityId_b;

        for (unsigned d=0; d < numDof_; ++d) {
          unsigned idx = numColEntities*numDof_+d;
          globalDofs_b[idx] = GID_(entityId_b, numDof_, d);
        }
        ++numColEntities;
    }
    globalDofs_b.resize(numColEntities*numDof_);

    // NOTE: 'Connections' should already include the self
    // pairings (where entity_a == entity_b) so we don't have
    // to worry about doing an insert on (globalRow_a, globalDofs_a),
    // etc.

    if (getDofStatus(entity_a) & ownedOrSharedMask) {
      //KOKKOS: small Loop noparallel insertGlobalIndices
      for (size_t d=0; d < numDof_; ++d) {
        const GlobalOrdinal globalRow_a = GID_(entityId_a, numDof_ , d);
        crsGraph.insertGlobalIndices(globalRow_a, globalDofs_b);
      }
    }

    for(unsigned j=0; j<numColEntities; ++j) {
      if (getDofStatus(entities_b[j]) & ownedOrSharedMask) {
        //KOKKOS: small Loop noparallel insertGlobalIndices
        for (size_t d=0; d < numDof_; ++d) {
          const GlobalOrdinal globalRow_b = GID_(entityIds_b[j], numDof_ , d);
          crsGraph.insertGlobalIndices(globalRow_b, globalDofs_a);
        }
      }
    }

    i += numColEntities;
    numColEntities = 0;
  }
}

void
TpetraLinearSystem::fill_entity_to_LID_mapping()
{
    const stk::mesh::BulkData& bulk = realm_.bulk_data();
    entityToLID_.assign(bulk.get_size_of_entity_index_space(), 200000000);
    const stk::mesh::BucketVector& nodeBuckets = bulk.buckets(stk::topology::NODE_RANK);
    for(const stk::mesh::Bucket* bptr : nodeBuckets) {
        const stk::mesh::Bucket& b = *bptr;
        const stk::mesh::EntityId* nodeIds = stk::mesh::field_data(*realm_.naluGlobalId_, b);
        for(size_t i=0; i<b.size(); ++i) {
            stk::mesh::Entity node = b[i];
            stk::mesh::EntityId nodeId = nodeIds[i];
            entityToLID_[node.local_offset()] = myLIDs_[nodeId];
        }
    }
}

void
TpetraLinearSystem::finalizeLinearSystem()
{
  ThrowRequire(inConstruction_);
  inConstruction_ = false;

  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  stk::mesh::MetaData & metaData = realm_.meta_data();

  const int this_mpi_rank = bulkData.parallel_rank();
  (void)this_mpi_rank;

  fill_entity_to_LID_mapping();

  ConnectionVec connectionVec;
  copy_kokkos_unordered_map_to_sorted_vector(connectionSetKK_, connectionVec);
  connectionSetKK_ = ConnectionSetKK(0);

  LinSys::RowLengths globallyOwnedRowLengths("rowLengths", globallyOwnedRowsMap_->getMyGlobalIndices().dimension(0));
  LinSys::RowLengths locallyOwnedRowLengths("rowLengths", ownedRowsMap_->getMyGlobalIndices().dimension(0));

  compute_graph_row_lengths(connectionVec, globallyOwnedRowLengths, locallyOwnedRowLengths);

  globallyOwnedGraph_ = Teuchos::rcp(new LinSys::Graph(globallyOwnedRowsMap_, ownedPlusGloballyOwnedRowsMap_, globallyOwnedRowLengths));
 
  insert_graph_connections(connectionVec, *globallyOwnedGraph_, DS_GloballyOwnedDOF);

  globallyOwnedGraph_->fillComplete();

  LinSys::Graph ownedPlusGloballyOwnedGraph(ownedRowsMap_, 8);
  ownedPlusGloballyOwnedGraph.doExport(*globallyOwnedGraph_, *exporter_, Tpetra::INSERT);
  ownedPlusGloballyOwnedGraph.fillComplete(ownedRowsMap_, ownedRowsMap_);

  // Add columns that are imported to the totalGids_ array
  const Teuchos::RCP<const LinSys::Map> & map = ownedPlusGloballyOwnedGraph.getColMap();
  for(size_t i=0; i<map->getNodeNumElements(); ++i) {
      const GlobalOrdinal gid = map->getGlobalElement(i);
      if(!ownedPlusGloballyOwnedRowsMap_->isNodeGlobalElement(gid))
        totalGids_.push_back(gid);
  }

  // This is the column map for the owned graph now
  const Teuchos::RCP<LinSys::Comm> tpetraComm = Teuchos::rcp(new LinSys::Comm(bulkData.parallel()));
  totalColsMap_ = Teuchos::rcp(new LinSys::Map(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), totalGids_, 1, tpetraComm, node_));
  ownedGraph_ = Teuchos::rcp(new LinSys::Graph(ownedRowsMap_, totalColsMap_, locallyOwnedRowLengths));

  insert_graph_connections(connectionVec, *ownedGraph_, DS_OwnedDOF);

  // add imported graph information
  {
    const LinSys::Map & rowMap = *ownedPlusGloballyOwnedGraph.getRowMap();
    const LinSys::Map & colMap = *ownedPlusGloballyOwnedGraph.getColMap();
    const size_t numRows = rowMap.getNodeNumElements();
    std::vector<GlobalOrdinal> newInd;
    for(size_t localRow=0; localRow<numRows; ++localRow) {
      const GlobalOrdinal row = rowMap.getGlobalElement(localRow);
      Teuchos::ArrayView<const LocalOrdinal> ind;
      ownedPlusGloballyOwnedGraph.getLocalRowView(localRow, ind);
      const size_t numInd = ind.size();
      newInd.resize(numInd);
      for(size_t j=0; j < numInd; ++j)
        {
          newInd[j] = colMap.getGlobalElement(ind[j]);
        }
      ownedGraph_->insertGlobalIndices(row, newInd);
    }
  }
  ownedGraph_->fillComplete(ownedRowsMap_, ownedRowsMap_);

  ownedMatrix_ = Teuchos::rcp(new LinSys::Matrix(ownedGraph_));
  globallyOwnedMatrix_ = Teuchos::rcp(new LinSys::Matrix(globallyOwnedGraph_));

  ownedLocalMatrix_ = ownedMatrix_->getLocalMatrix();
  globallyOwnedLocalMatrix_ = globallyOwnedMatrix_->getLocalMatrix();

  ownedRhs_ = Teuchos::rcp(new LinSys::Vector(ownedRowsMap_));
  globallyOwnedRhs_ = Teuchos::rcp(new LinSys::Vector(globallyOwnedRowsMap_));

  ownedLocalRhs_ = ownedRhs_->getLocalView<sierra::nalu::HostSpace>();
  globallyOwnedLocalRhs_ = globallyOwnedRhs_->getLocalView<sierra::nalu::HostSpace>();

  sln_ = Teuchos::rcp(new LinSys::Vector(ownedRowsMap_));

  const int nDim = metaData.spatial_dimension();

  Teuchos::RCP<LinSys::MultiVector> coords 
    = Teuchos::RCP<LinSys::MultiVector>(new LinSys::MultiVector(sln_->getMap(), nDim));

  TpetraLinearSolver *linearSolver = reinterpret_cast<TpetraLinearSolver *>(linearSolver_);

  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  if (linearSolver->activeMueLu())
    copy_stk_to_tpetra(coordinates, coords);

  linearSolver->setupLinearSolver(sln_, ownedMatrix_, ownedRhs_, coords);
}

void
TpetraLinearSystem::zeroSystem()
{
  ThrowRequire(!ownedMatrix_.is_null());
  ThrowRequire(!globallyOwnedMatrix_.is_null());
  ThrowRequire(!globallyOwnedRhs_.is_null());
  ThrowRequire(!ownedRhs_.is_null());

  globallyOwnedMatrix_->resumeFill();
  ownedMatrix_->resumeFill();

  globallyOwnedMatrix_->setAllToScalar(0);
  ownedMatrix_->setAllToScalar(0);
  globallyOwnedRhs_->putScalar(0);
  ownedRhs_->putScalar(0);

  sln_->putScalar(0);
}

namespace
{
  template <typename RowViewType>
  void sum_into_row (
    RowViewType row_view,
    int numCols,
    const int* localIds,
    const int* sort_permutation,
    const double* input_values)
  {
    constexpr bool forceAtomic = !std::is_same<sierra::nalu::DeviceSpace, Kokkos::Serial>::value;
    const LocalOrdinal length = row_view.length;

    LocalOrdinal offset = 0;
    for (int j = 0; j < numCols; ++j) {
      const LocalOrdinal perm_index = sort_permutation[j];
      const LocalOrdinal cur_local_column_idx = localIds[j];

      // since the columns are sorted, we pass through the column idxs once,
      // updating the offset as we go
      while (row_view.colidx(offset) != cur_local_column_idx && offset < length) {
        ++offset;
      }

      if (offset < length) {
        ThrowAssertMsg(std::isfinite(input_values[perm_index]), "Inf or NAN lhs");
        if (forceAtomic) {
          Kokkos::atomic_add(&(row_view.value(offset)), input_values[perm_index]);
        }
        else {
          row_view.value(offset) += input_values[perm_index];
        }
      }
    }
  }
}

void
TpetraLinearSystem::sumInto(
      unsigned numEntities,
      const stk::mesh::Entity* entities,
      const SharedMemView<const double*> & rhs,
      const SharedMemView<const double**> & lhs,
      const SharedMemView<int*> & localIds,
      const SharedMemView<int*> & sortPermutation,
      const char * trace_tag)
{
  constexpr bool forceAtomic = !std::is_same<sierra::nalu::DeviceSpace, Kokkos::Serial>::value;

  ThrowAssertMsg(lhs.is_contiguous(), "LHS assumed contiguous");
  ThrowAssertMsg(rhs.is_contiguous(), "RHS assumed contiguous");
  ThrowAssertMsg(localIds.is_contiguous(), "localIds assumed contiguous");
  ThrowAssertMsg(sortPermutation.is_contiguous(), "sortPermutation assumed contiguous");

  const int n_obj = numEntities;
  const int numRows = n_obj * numDof_;

  for(int i = 0; i < n_obj; i++) {
    const stk::mesh::Entity entity = entities[i];
    const LocalOrdinal localOffset = entityToLID_[entity.local_offset()] * numDof_;
    for(size_t d=0; d < numDof_; ++d) {
      size_t lid = i*numDof_ + d;
      localIds[lid] = localOffset + d;
    }
  }

  for (int i = 0; i < numRows; ++i) {
    sortPermutation[i] = i;
  }
  Tpetra::Details::shellSortKeysAndValues(localIds.data(), sortPermutation.data(), numRows);

  for (int r = 0; r < numRows; r++) {
    const LocalOrdinal localId = localIds[r];
    const LocalOrdinal cur_perm_index = sortPermutation[r];
    const double* const cur_lhs = &lhs(cur_perm_index, 0);
    const double cur_rhs = rhs[cur_perm_index];
    ThrowAssertMsg(std::isfinite(cur_rhs), "Inf or NAN rhs");

    if(localId < maxOwnedRowId_) {
      sum_into_row(ownedLocalMatrix_.row(localId), numRows, localIds.data(), sortPermutation.data(), cur_lhs);
      if (forceAtomic) {
        Kokkos::atomic_add(&ownedLocalRhs_(localId,0), cur_rhs);
      }
      else {
        ownedLocalRhs_(localId,0) += cur_rhs;
      }
    }
    else if (localId < maxGloballyOwnedRowId_) {
      const LocalOrdinal actualLocalId = localId - maxOwnedRowId_;
      sum_into_row(globallyOwnedLocalMatrix_.row(actualLocalId), numRows,
        localIds.data(), sortPermutation.data(), cur_lhs);

      if (forceAtomic) {
        Kokkos::atomic_add(&globallyOwnedLocalRhs_(actualLocalId,0), cur_rhs);
      }
      else {
        globallyOwnedLocalRhs_(actualLocalId,0) += cur_rhs;
      }
    }
  }
}

void
TpetraLinearSystem::sumInto(
  const std::vector<stk::mesh::Entity> & entities,
  std::vector<int> &scratchIds,
  std::vector<double> &scratchVals,
  const std::vector<double> & rhs,
  const std::vector<double> & lhs,
  const char *trace_tag
  )
{
  const size_t n_obj = entities.size();
  const unsigned numRows = n_obj * numDof_;

  ThrowAssert(numRows == rhs.size());
  ThrowAssert(numRows*numRows == lhs.size());

  scratchIds.resize(numRows);
  sortPermutation_.resize(numRows);
  for(size_t i = 0; i < n_obj; i++) {
    const stk::mesh::Entity entity = entities[i];
    const LocalOrdinal localOffset = entityToLID_[entity.local_offset()] * numDof_;
    for(size_t d=0; d < numDof_; ++d) {
      size_t lid = i*numDof_ + d;
      scratchIds[lid] = localOffset + d;
    }
  }

  for (unsigned i = 0; i < numRows; ++i) {
    sortPermutation_[i] = i;
  }
  Tpetra::Details::shellSortKeysAndValues(scratchIds.data(), sortPermutation_.data(), (int)numRows);

  for (unsigned r = 0; r < numRows; r++) {
    const LocalOrdinal localId = scratchIds[r];
    const LocalOrdinal cur_perm_index = sortPermutation_[r];
    const double* const cur_lhs = &lhs[cur_perm_index*numRows];
    const double cur_rhs = rhs[cur_perm_index];
    ThrowAssertMsg(std::isfinite(cur_rhs), "Invalid rhs");

    if(localId < maxOwnedRowId_) {
      sum_into_row(ownedLocalMatrix_.row(localId), numRows, scratchIds.data(), sortPermutation_.data(), cur_lhs);
      ownedLocalRhs_(localId,0) += cur_rhs;
    }
    else if (localId < maxGloballyOwnedRowId_) {
      const LocalOrdinal actualLocalId = localId - maxOwnedRowId_;
      sum_into_row(globallyOwnedLocalMatrix_.row(actualLocalId), numRows,
        scratchIds.data(), sortPermutation_.data(), cur_lhs);

      globallyOwnedLocalRhs_(actualLocalId,0) += cur_rhs;
    }
  }
}

void
TpetraLinearSystem::applyDirichletBCs(
  stk::mesh::FieldBase * solutionField,
  stk::mesh::FieldBase * bcValuesField,
  const stk::mesh::PartVector & parts,
  const unsigned beginPos,
  const unsigned endPos)
{
  stk::mesh::MetaData & metaData = realm_.meta_data();

  double adbc_time = -NaluEnv::self().nalu_time();

  const stk::mesh::Selector selector 
    = (metaData.locally_owned_part() | metaData.globally_shared_part())
    & stk::mesh::selectUnion(parts)
    & stk::mesh::selectField(*solutionField) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, selector );

  const bool internalMatrixIsSorted = true;
  int nbc=0;
  for(const stk::mesh::Bucket* bptr : buckets) {
    const stk::mesh::Bucket & b = *bptr;

    const unsigned fieldSize = field_bytes_per_entity(*solutionField, b) / sizeof(double);
    ThrowRequire(fieldSize == numDof_);

    const stk::mesh::Bucket::size_type length   = b.size();
    const double * solution = (double*)stk::mesh::field_data(*solutionField, *b.begin());
    const double * bcValues = (double*)stk::mesh::field_data(*bcValuesField, *b.begin());

    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const double> values;
    std::vector<double> new_values;

    for (stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const stk::mesh::Entity entity = b[k];
      const stk::mesh::EntityId naluId = *stk::mesh::field_data(*realm_.naluGlobalId_, entity);
      const LocalOrdinal localIdOffset = lookup_myLID(myLIDs_, naluId, "applyDirichletBCs") * numDof_;

      for(unsigned d=beginPos; d < endPos; ++d) {
        const LocalOrdinal localId = localIdOffset + d;
        const bool useOwned = localId < maxOwnedRowId_;
        const LocalOrdinal actualLocalId = useOwned ? localId : localId - maxOwnedRowId_;
        Teuchos::RCP<LinSys::Matrix> matrix = useOwned ? ownedMatrix_ : globallyOwnedMatrix_;
        const LinSys::Matrix::local_matrix_type& local_matrix = useOwned ? ownedLocalMatrix_ : globallyOwnedLocalMatrix_;

        if(localId > maxGloballyOwnedRowId_) {
          std::cerr << "localId > maxGloballyOwnedRowId_:: localId= " << localId << " maxGloballyOwnedRowId_= " << maxGloballyOwnedRowId_ << " useOwned = " << (localId < maxOwnedRowId_ ) << std::endl;
          throw std::runtime_error("logic error: localId > maxGloballyOwnedRowId_");
        }

        // Adjust the LHS

        const double diagonal_value = useOwned ? 1.0 : 0.0;

        matrix->getLocalRowView(actualLocalId, indices, values);
        const size_t rowLength = values.size();
        if (rowLength > 0) {
          new_values.resize(rowLength);
          for(size_t i=0; i < rowLength; ++i) {
              new_values[i] = (indices[i] == localId) ? diagonal_value : 0;
          }
          local_matrix.replaceValues(actualLocalId, &indices[0], rowLength, new_values.data(), internalMatrixIsSorted);
        }

        // Replace the RHS residual with (desired - actual)
        Teuchos::RCP<LinSys::Vector> rhs = useOwned ? ownedRhs_: globallyOwnedRhs_;
        const double bc_residual = useOwned ? (bcValues[k*fieldSize + d] - solution[k*fieldSize + d]) : 0.0;
        rhs->replaceLocalValue(actualLocalId, bc_residual);
        ++nbc;
      }
    }
  }
  adbc_time += NaluEnv::self().nalu_time();
}

void
TpetraLinearSystem::prepareConstraints(
  const unsigned beginPos,
  const unsigned endPos)
{
  Teuchos::ArrayView<const LocalOrdinal> indices;
  Teuchos::ArrayView<const double> values;
  std::vector<double> new_values;

  const bool internalMatrixIsSorted = true;

  // iterate oversetInfoVec_
  std::vector<OversetInfo *>::iterator ii;
  //KOKKOS: Loop noparallel RCP Vector Matrix replaceValues
  for( ii=realm_.oversetManager_->oversetInfoVec_.begin();
       ii!=realm_.oversetManager_->oversetInfoVec_.end(); ++ii ) {

    // extract orphan node and global id; process both owned and shared
    stk::mesh::Entity orphanNode = (*ii)->orphanNode_;
    const stk::mesh::EntityId naluId = *stk::mesh::field_data(*realm_.naluGlobalId_, orphanNode);
    const LocalOrdinal localIdOffset = lookup_myLID(myLIDs_, naluId, "prepareConstraints") * numDof_;

    //KOKKOS: Nested Loop noparallel RCP Vector Matrix replaceValues
    for(unsigned d=beginPos; d < endPos; ++d) {
      const LocalOrdinal localId = localIdOffset + d;
      const bool useOwned = localId < maxOwnedRowId_;
      const LocalOrdinal actualLocalId = useOwned ? localId : localId - maxOwnedRowId_;
      Teuchos::RCP<LinSys::Matrix> matrix = useOwned ? ownedMatrix_ : globallyOwnedMatrix_;
      const LinSys::Matrix::local_matrix_type& local_matrix = matrix->getLocalMatrix();
      
      if ( localId > maxGloballyOwnedRowId_) {
        throw std::runtime_error("logic error: localId > maxGloballyOwnedRowId_");
      }
      
      // Adjust the LHS; full row is perfectly zero
      matrix->getLocalRowView(actualLocalId, indices, values);
      const size_t rowLength = values.size();
      if (rowLength > 0) {
        new_values.resize(rowLength);
        for(size_t i=0; i < rowLength; ++i) {
          new_values[i] = 0.0;
        }
        local_matrix.replaceValues(actualLocalId, &indices[0], rowLength, new_values.data(), internalMatrixIsSorted);
      }
      
      // Replace the RHS residual with zero
      Teuchos::RCP<LinSys::Vector> rhs = useOwned ? ownedRhs_: globallyOwnedRhs_;
      const double bc_residual = 0.0;
      rhs->replaceLocalValue(actualLocalId, bc_residual);
    }
  }
}

void
TpetraLinearSystem::resetRows(
  const std::vector<stk::mesh::Entity> nodeList,
  const unsigned beginPos,
  const unsigned endPos)
{
  Teuchos::ArrayView<const LocalOrdinal> indices;
  Teuchos::ArrayView<const double> values;
  std::vector<double> new_values;
  constexpr double rhs_residual = 0.0;
  const bool internalMatrixIsSorted = true;

  for (auto node: nodeList) {
    const auto naluId = *stk::mesh::field_data(*realm_.naluGlobalId_, node);
    const LocalOrdinal localIdOffset = lookup_myLID(
      myLIDs_, naluId, "resetRows") * numDof_;

    for (unsigned d=beginPos; d < endPos; ++d) {
      const LocalOrdinal localId = localIdOffset + d;
      const bool useOwned = (localId < maxOwnedRowId_);
      const LocalOrdinal actualLocalId =
        useOwned ? localId : (localId - maxOwnedRowId_);
      Teuchos::RCP<LinSys::Matrix> matrix =
        useOwned ? ownedMatrix_ : globallyOwnedMatrix_;
      const LinSys::Matrix::local_matrix_type& local_matrix = matrix->getLocalMatrix();

      if (localId > maxGloballyOwnedRowId_) {
        throw std::runtime_error("logic error: localId > maxGloballyOwnedRowId");
      }

      // Adjust the LHS; zero out all entries (including diagonal)
      matrix->getLocalRowView(actualLocalId, indices, values);
      const size_t rowLength = values.size();
      if (rowLength > 0) {
        new_values.resize(rowLength);
        for (size_t i=0; i < rowLength; i++) {
          new_values[i] = 0.0;
        }
        local_matrix.replaceValues(actualLocalId, &indices[0], rowLength, new_values.data(), internalMatrixIsSorted);
      }

      // Replace RHS residual entry = 0.0
      Teuchos::RCP<LinSys::Vector> rhs =
        useOwned ? ownedRhs_ : globallyOwnedRhs_;
      rhs->replaceLocalValue(actualLocalId, rhs_residual);
    }
  }
}

void
TpetraLinearSystem::loadComplete()
{
  // LHS
  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::parameterList ();
  params->set("No Nonlocal Changes", true);
  bool do_params=false;

  if (do_params)
    globallyOwnedMatrix_->fillComplete(params);
  else
    globallyOwnedMatrix_->fillComplete();

  ownedMatrix_->doExport(*globallyOwnedMatrix_, *exporter_, Tpetra::ADD);
  if (do_params)
    ownedMatrix_->fillComplete(params);
  else
    ownedMatrix_->fillComplete();

  // RHS
  ownedRhs_->doExport(*globallyOwnedRhs_, *exporter_, Tpetra::ADD);
}

int
TpetraLinearSystem::solve(
  stk::mesh::FieldBase * linearSolutionField)
{

  TpetraLinearSolver *linearSolver = reinterpret_cast<TpetraLinearSolver *>(linearSolver_);

  if ( realm_.debug() ) {
    checkForNaN(true);
    if (checkForZeroRow(true, false, true)) {
      throw std::runtime_error("ERROR checkForZeroRow in solve()");
    }
  }
   
  if (linearSolver->getConfig()->getWriteMatrixFiles()) {
    writeToFile(eqSysName_.c_str());
    writeToFile(eqSysName_.c_str(), false);
  }

  double solve_time = -NaluEnv::self().nalu_time();

  int iters;
  double finalResidNorm;
  
  // memory diagnostic
  if ( realm_.get_activate_memory_diagnostic() ) {
    NaluEnv::self().naluOutputP0() << "NaluMemory::TpetraLinearSystem::solve() PreSolve: " << eqSysName_ << std::endl;
    realm_.provide_memory_summary();
  }

  const int status = linearSolver->solve(
      sln_,
      iters,
      finalResidNorm,
      realm_.isFinalOuterIter_);

  solve_time += NaluEnv::self().nalu_time();

  if (linearSolver->getConfig()->getWriteMatrixFiles()) {
    writeSolutionToFile(eqSysName_.c_str());
    ++writeCounter_;
  }

  copy_tpetra_to_stk(sln_, linearSolutionField);
  sync_field(linearSolutionField);

  // computeL2 norm
  const double norm2 = ownedRhs_->norm2();

  // save off solver info
  linearSolveIterations_ = iters;
  nonLinearResidual_ = realm_.l2Scaling_*norm2;
  linearResidual_ = finalResidNorm;
   
  if ( eqSys_->firstTimeStepSolve_ )
    firstNonLinearResidual_ = nonLinearResidual_;
  scaledNonLinearResidual_ = nonLinearResidual_/std::max(std::numeric_limits<double>::epsilon(), firstNonLinearResidual_);

  if ( provideOutput_ ) {
    const int nameOffset = eqSysName_.length()+8;
    NaluEnv::self().naluOutputP0()
      << std::setw(nameOffset) << std::right << eqSysName_
      << std::setw(32-nameOffset)  << std::right << iters
      << std::setw(18) << std::right << finalResidNorm
      << std::setw(15) << std::right << nonLinearResidual_
      << std::setw(14) << std::right << scaledNonLinearResidual_ << std::endl;
  }

  eqSys_->firstTimeStepSolve_ = false;

  return status;
}

void
TpetraLinearSystem::checkForNaN(bool useOwned)
{
  Teuchos::RCP<LinSys::Matrix> matrix = useOwned ? ownedMatrix_ : globallyOwnedMatrix_;
  Teuchos::RCP<LinSys::Vector> rhs = useOwned ? ownedRhs_ : globallyOwnedRhs_;

  Teuchos::ArrayView<const LocalOrdinal> indices;
  Teuchos::ArrayView<const double> values;

  size_t n = matrix->getRowMap()->getNodeNumElements();
  for(size_t i=0; i<n; ++i) {

    matrix->getLocalRowView(i, indices, values);
    const size_t rowLength = values.size();
    for(size_t k=0; k < rowLength; ++k) {
      if (values[k] != values[k])	{
        std::cerr << "LHS NaN: " << i << std::endl;
        throw std::runtime_error("bad LHS");
      }
    }
  }

  Teuchos::ArrayRCP<const Scalar> rhs_data = rhs->getData();
  n = rhs_data.size();
  for(size_t i=0; i<n; ++i) {
    if (rhs_data[i] != rhs_data[i]) {
      std::cerr << "rhs NaN: " << i << std::endl;
      throw std::runtime_error("bad rhs");
    }
  }
}

bool
TpetraLinearSystem::checkForZeroRow(bool useOwned, bool doThrow, bool doPrint)
{
  Teuchos::RCP<LinSys::Matrix> matrix = useOwned ? ownedMatrix_ : globallyOwnedMatrix_;
  Teuchos::RCP<LinSys::Vector> rhs = useOwned ? ownedRhs_ : globallyOwnedRhs_;
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  Teuchos::ArrayView<const LocalOrdinal> indices;
  Teuchos::ArrayView<const double> values;

  size_t nrowG = matrix->getRangeMap()->getGlobalNumElements();
  size_t n = matrix->getRowMap()->getNodeNumElements();
  GlobalOrdinal max_gid = 0, g_max_gid=0;
  //KOKKOS: Loop parallel reduce
  kokkos_parallel_for("Nalu::TpetraLinearSystem::checkForZeroRowA", n, [&] (const size_t& i) {
    GlobalOrdinal gid = matrix->getGraph()->getRowMap()->getGlobalElement(i);
    max_gid = std::max(gid, max_gid);
  });
  stk::all_reduce_max(bulkData.parallel(), &max_gid, &g_max_gid, 1);

  nrowG = g_max_gid+1;
  std::vector<double> local_row_sums(nrowG, 0.0);
  std::vector<int> local_row_exists(nrowG, 0);
  std::vector<double> global_row_sums(nrowG, 0.0);
  std::vector<int> global_row_exists(nrowG, 0);

  for(size_t i=0; i<n; ++i) {
    GlobalOrdinal gid = matrix->getGraph()->getRowMap()->getGlobalElement(i);
    matrix->getLocalRowView(i, indices, values);
    const size_t rowLength = values.size();
    double row_sum = 0.0;
    for(size_t k=0; k < rowLength; ++k) {
      row_sum += std::abs(values[k]);
    }
    if (gid-1 >= (GlobalOrdinal)local_row_sums.size() || gid <= 0) {
      std::cerr << "gid= " << gid << " nrowG= " << nrowG << std::endl;
      throw std::runtime_error("bad gid");
    }
    local_row_sums[gid-1] = row_sum;
    local_row_exists[gid-1] = 1;
  }

  stk::all_reduce_sum(bulkData.parallel(), &local_row_sums[0], &global_row_sums[0], (unsigned)nrowG);
  stk::all_reduce_max(bulkData.parallel(), &local_row_exists[0], &global_row_exists[0], (unsigned)nrowG);

  bool found=false;
  //KOKKOS: Loop parallel
  kokkos_parallel_for("Nalu::TpetraLinearSystem::checkForZeroRowC", nrowG, [&] (const size_t& ii) {
    double row_sum = global_row_sums[ii];
    if (global_row_exists[ii] && bulkData.parallel_rank() == 0 && row_sum < 1.e-10) {
      found = true;
      GlobalOrdinal gid = ii+1;
      stk::mesh::EntityId nid = GLOBAL_ENTITY_ID(gid, numDof_);
      stk::mesh::Entity node = bulkData.get_entity(stk::topology::NODE_RANK, nid);
      stk::mesh::EntityId naluGlobalId;
      if (bulkData.is_valid(node)) naluGlobalId = *stk::mesh::field_data(*realm_.naluGlobalId_, node);

      int idof = GLOBAL_ENTITY_ID_IDOF(gid, numDof_);
      GlobalOrdinal GID_check = GID_(nid, numDof_, idof);
      if (doPrint) {

        double dualVolume = -1.0;

        std::cout << "P[" << bulkData.parallel_rank() << "] LHS zero: " << ii
                  << " GID= " << gid << " GID_check= " << GID_check << " nid= " << nid
                  << " naluGlobalId " << naluGlobalId << " is_valid= " << bulkData.is_valid(node)
                  << " idof= " << idof << " numDof_= " << numDof_
                  << " row_sum= " << row_sum
                  << " dualVolume= " << dualVolume
                  << std::endl;
        NaluEnv::self().naluOutputP0() << "P[" << bulkData.parallel_rank() << "] LHS zero: " << ii
                        << " GID= " << gid << " GID_check= " << GID_check << " nid= " << nid
                        << " naluGlobalId " << naluGlobalId << " is_valid= " << bulkData.is_valid(node)
                        << " idof= " << idof << " numDof_= " << numDof_
                        << " row_sum= " << row_sum
                        << " dualVolume= " << dualVolume
                        << std::endl;
      }
    }
  });

  if (found && doThrow) {
    throw std::runtime_error("bad zero row LHS");
  }
  return found;
}

void
TpetraLinearSystem::writeToFile(const char * base_filename, bool useOwned)
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const unsigned p_rank = bulkData.parallel_rank();
  const unsigned p_size = bulkData.parallel_size();

  Teuchos::RCP<LinSys::Matrix> matrix = useOwned ? ownedMatrix_ : globallyOwnedMatrix_;
  Teuchos::RCP<LinSys::Vector> rhs = useOwned ? ownedRhs_ : globallyOwnedRhs_;

  const int currentCount = writeCounter_;

  if (1)
    {
      std::ostringstream osLhs;
      std::ostringstream osRhs;
      osLhs << base_filename << "-" << (useOwned ? "O-":"G-") << currentCount << ".mm." << p_size; // A little hacky but whatever
      osRhs << base_filename << "-" << (useOwned ? "O-":"G-") << currentCount << ".rhs." << p_size; // A little hacky but whatever

      Tpetra::MatrixMarket::Writer<LinSys::Matrix>::writeSparseFile(osLhs.str().c_str(), matrix,
                                                                    eqSysName_, std::string("Tpetra matrix for: ")+eqSysName_, true);
      typedef Tpetra::MatrixMarket::Writer<LinSys::Matrix> writer_type;
      if (useOwned) writer_type::writeDenseFile (osRhs.str().c_str(), rhs);
    }

  if (1)
    {
      std::ostringstream osLhs;
      std::ostringstream osGra;
      std::ostringstream osRhs;

      osLhs << base_filename << "-" << (useOwned ? "O-":"G-") << currentCount << ".mm." << p_size << "." << p_rank; // A little hacky but whatever
      osGra << base_filename << "-" << (useOwned ? "O-":"G-") << currentCount << ".gra." << p_size << "." << p_rank; // A little hacky but whatever
      osRhs << base_filename << "-" << (useOwned ? "O-":"G-") << currentCount << ".rhs." << p_size << "." << p_rank; // A little hacky but whatever

      //Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
#define DUMP(A)  do {                                                   \
        out << "\n\n===============================================================================================\n"; \
        out << "===============================================================================================\n"; \
        out << "P[" << p_rank << "] writeToFile:: " #A "= " << "\n---------------------------\n" ; \
        out << Teuchos::describe(*A,Teuchos::VERB_EXTREME) << "\n";     \
        out << "===============================================================================================\n"; \
        out << "===============================================================================================\n\n\n"; \
      } while(0)

      {
        std::ostringstream out;
        DUMP(matrix);
        std::ofstream fout;
        fout.open (osLhs.str().c_str());
        fout << out.str() << std::endl;
      }

      {
        std::ostringstream out;
        DUMP(matrix->getGraph());
        std::ofstream fout;
        fout.open (osGra.str().c_str());
        fout << out.str() << std::endl;
      }

      {
        std::ostringstream out;
        DUMP(rhs);
        std::ofstream fout;
        fout.open (osRhs.str().c_str());
        fout << out.str() << std::endl;
      }


#undef DUMP

    }

}

void
TpetraLinearSystem::printInfo(bool useOwned)
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const unsigned p_rank = bulkData.parallel_rank();

  Teuchos::RCP<LinSys::Matrix> matrix = useOwned ? ownedMatrix_ : globallyOwnedMatrix_;
  Teuchos::RCP<LinSys::Vector> rhs = useOwned ? ownedRhs_ : globallyOwnedRhs_;

  if (p_rank == 0) {
    std::cout << "\nMatrix for EqSystem: " << eqSysName_ << " :: N N NZ= " << matrix->getRangeMap()->getGlobalNumElements()
              << " "
              << matrix->getDomainMap()->getGlobalNumElements()
              << " "
              << matrix->getGlobalNumEntries()
              << std::endl;
    NaluEnv::self().naluOutputP0() << "\nMatrix for system: " << eqSysName_ << " :: N N NZ= " << matrix->getRangeMap()->getGlobalNumElements()
                                   << " "
                                   << matrix->getDomainMap()->getGlobalNumElements()
                                   << " "
                                   << matrix->getGlobalNumEntries()
                                   << std::endl;
  }
}

void
TpetraLinearSystem::writeSolutionToFile(const char * base_filename, bool useOwned)
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const unsigned p_rank = bulkData.parallel_rank();
  const unsigned p_size = bulkData.parallel_size();

  Teuchos::RCP<LinSys::Vector> sln = sln_;
  const int currentCount = writeCounter_;

  if (1)
    {
      std::ostringstream osSln;
      osSln << base_filename << "-" << (useOwned ? "O-":"G-") << currentCount << ".sln." << p_size; // A little hacky but whatever

      typedef Tpetra::MatrixMarket::Writer<LinSys::Matrix> writer_type;
      if (useOwned) writer_type::writeDenseFile (osSln.str().c_str(), sln);
    }

  if (1)
    {
      std::ostringstream osSln;

      osSln << base_filename << "-" << "O-" << currentCount << ".sln." << p_size << "." << p_rank; // A little hacky but whatever

#define DUMP(A)  do {                                                   \
        out << "\n\n===============================================================================================\n"; \
        out << "===============================================================================================\n"; \
        out << "P[" << p_rank << "] writeToFile:: " #A "= " << "\n---------------------------\n" ; \
        out << Teuchos::describe(*A,Teuchos::VERB_EXTREME) << "\n";     \
        out << "===============================================================================================\n"; \
        out << "===============================================================================================\n\n\n"; \
      } while(0)

      {
        std::ostringstream out;
        DUMP(sln);
        std::ofstream fout;
        fout.open (osSln.str().c_str());
        fout << out.str() << std::endl;
      }


#undef DUMP

    }

}

void
TpetraLinearSystem::copy_tpetra_to_stk(
  const Teuchos::RCP<LinSys::Vector> tpetraField,
  stk::mesh::FieldBase * stkField)
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  stk::mesh::MetaData & metaData = realm_.meta_data();

  ThrowAssert(!tpetraField.is_null());
  ThrowAssert(stkField);
  const LinSys::ConstOneDVector & tpetraVector = tpetraField->get1dView();

  const unsigned p_rank = bulkData.parallel_rank();

  const stk::mesh::Selector selector = stk::mesh::selectField(*stkField)
    & metaData.locally_owned_part() 
    & !(stk::mesh::selectUnion(realm_.get_slave_part_vector()))
    & !(realm_.get_inactive_selector());
  
  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets(stk::topology::NODE_RANK, selector);

  for (size_t ib=0; ib < buckets.size(); ++ib) {
    stk::mesh::Bucket & b = *buckets[ib];

    const unsigned fieldSize = field_bytes_per_entity(*stkField, b) / sizeof(double);
    ThrowRequire(fieldSize == numDof_);

    const stk::mesh::Bucket::size_type length = b.size();
    double * stkFieldPtr = (double*)stk::mesh::field_data(*stkField, *b.begin());
    const stk::mesh::EntityId *naluGlobalId = stk::mesh::field_data(*realm_.naluGlobalId_, *b.begin());
    for (stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      stk::mesh::Entity node = b[k];
      const LocalOrdinal localIdOffset = entityToLID_[node.local_offset()] * numDof_;
      for(unsigned d=0; d < fieldSize; ++d) {
        const LocalOrdinal localId = localIdOffset + d;
        bool useOwned = true;
        LocalOrdinal actualLocalId = localId;
        if(localId >= maxOwnedRowId_) {
          actualLocalId = localId - maxOwnedRowId_;
          useOwned = false;
        }

        if (!useOwned) {
          stk::mesh::EntityId naluId = naluGlobalId[k];
          stk::mesh::EntityId stkId = bulkData.identifier(node);
          std::cout << "P[" << p_rank << "] useOwned = " << useOwned << " localId = " << localId << " maxOwnedRowId_= " << maxOwnedRowId_ << " actualLocalId= " << actualLocalId
                    << " naluGlobalId= " << naluGlobalId[k] << " stkId= " << stkId << " naluId= " << naluId << std::endl;
        }
        ThrowRequire(useOwned);

        const size_t stkIndex = k*numDof_ + d;
        if (useOwned){
          stkFieldPtr[stkIndex] = tpetraVector[localId];
        }
      }
    }
  }
}

int getDofStatus_impl(stk::mesh::Entity node, const Realm& realm)
{
  const stk::mesh::BulkData & bulkData = realm.bulk_data();

  const stk::mesh::Bucket & b = bulkData.bucket(node);
  const bool entityIsOwned = b.owned();
  const bool entityIsShared = b.shared();
  const bool entityIsGhosted = !entityIsOwned && !entityIsShared;

  bool has_non_matching_boundary_face_alg = realm.has_non_matching_boundary_face_alg();
  bool hasPeriodic = realm.hasPeriodic_;

  if (realm.hasPeriodic_ && realm.has_non_matching_boundary_face_alg()) {
    has_non_matching_boundary_face_alg = false;
    hasPeriodic = false;

    stk::mesh::Selector perSel = stk::mesh::selectUnion(realm.allPeriodicInteractingParts_);
    stk::mesh::Selector nonConfSel = stk::mesh::selectUnion(realm.allNonConformalInteractingParts_);
    //std::cout << "nonConfSel= " << nonConfSel << std::endl;

    for (auto part : b.supersets()) {
      if (perSel(*part)) {
        hasPeriodic = true;
      }
      if (nonConfSel(*part)) {
        has_non_matching_boundary_face_alg = true;
      }
    }
  }

  //std::cerr << "has_non_matching_boundary_face_alg= " << has_non_matching_boundary_face_alg << " hasPeriodic= " << hasPeriodic << std::endl;

  if (has_non_matching_boundary_face_alg && hasPeriodic) {
    std::ostringstream ostr;
    ostr << "node id= " << realm.bulkData_->identifier(node);
    throw std::logic_error("not ready for primetime to combine periodic and non-matching algorithm on same node: "+ostr.str());
  }

  // simple case
  if (!hasPeriodic && !has_non_matching_boundary_face_alg) {
    if (entityIsGhosted)
      return DS_GhostedDOF;
    if (entityIsOwned)
      return DS_OwnedDOF;
    if (entityIsShared && !entityIsOwned)
      return DS_GloballyOwnedDOF;
  }

  if (has_non_matching_boundary_face_alg) {
    if (entityIsOwned)
      return DS_OwnedDOF;
    //if (entityIsShared && !entityIsOwned) {
    if (!entityIsOwned && (entityIsGhosted || entityIsShared)){
      return DS_GloballyOwnedDOF;
    }
    // maybe return DS_GhostedDOF if entityIsGhosted
  }

  if (hasPeriodic) {
    const stk::mesh::EntityId stkId = bulkData.identifier(node);
    const stk::mesh::EntityId naluId = *stk::mesh::field_data(*realm.naluGlobalId_, node);

    // bool for type of ownership for this node
    const bool nodeOwned = bulkData.bucket(node).owned();
    const bool nodeShared = bulkData.bucket(node).shared();
    const bool nodeGhosted = !nodeOwned && !nodeShared;

    // really simple here.. ghosted nodes never part of the matrix
    if ( nodeGhosted ) {
      return DS_GhostedDOF;
    }

    // bool to see if this is possibly a periodic node
    const bool isSlaveNode = (stkId != naluId);

    if (!isSlaveNode) {
      if (nodeOwned)
        return DS_OwnedDOF;
      else if( nodeShared )
        return DS_GloballyOwnedDOF;
      else
        return DS_GhostedDOF;
    }
    else {
      // I am a slave node.... get the master entity
      stk::mesh::Entity masterEntity = bulkData.get_entity(stk::topology::NODE_RANK, naluId);
      if ( bulkData.is_valid(masterEntity)) {
        const bool masterEntityOwned = bulkData.bucket(masterEntity).owned();
        const bool masterEntityShared = bulkData.bucket(masterEntity).shared();
        if (masterEntityOwned)
          return DS_SkippedDOF | DS_OwnedDOF;
        if (masterEntityShared)
          return DS_SkippedDOF | DS_GloballyOwnedDOF;
        else
          return DS_GloballyOwnedDOF;
      }
      else {
        return DS_GloballyOwnedDOF;
      }
    }
  }

  // still got here? problem...
  if (1)
    throw std::logic_error("bad status2");

  return DS_SkippedDOF;
}

} // namespace nalu
} // namespace Sierra
