/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <TpetraLinearSystem.h>
#include <ContactInfo.h>
#include <ContactManager.h>
#include <NonConformalInfo.h>
#include <NonConformalManager.h>
#include <FieldTypeDef.h>
#include <HaloInfo.h>
#include <DgInfo.h>
#include <Realm.h>
#include <PeriodicManager.h>
#include <Simulation.h>
#include <LinearSolver.h>
#include <master_element/MasterElement.h>
#include <NaluEnv.h>

// overset
#include <overset/OversetManager.h>
#include <overset/OversetInfo.h>

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

#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Tpetra_MatrixIO.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <set>
#include <limits>

#include <sstream>

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
  const std::string & name,
  LinearSolver * linearSolver)
  : LinearSystem(realm, numDof, name, linearSolver)
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
};

size_t TpetraLinearSystem::lookup_myLID(MyLIDMapType& myLIDs, stk::mesh::EntityId entityId, const std::string& msg, stk::mesh::Entity entity)
{
  return myLIDs[entityId];
}

// determines whether the node is to be put into which map/graph/matrix
// FIXME - note that the DOFStatus enum can be Or'd together if need be to
//   distinguish ever more complicated situations, for example, a DOF that
//   is both owned and ghosted: OwnedDOF | GhostedDOF
int TpetraLinearSystem::getDofStatus(stk::mesh::Entity node)
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const unsigned p_rank = bulkData.parallel_rank();
  (void)p_rank;

  const stk::mesh::Bucket & b = bulkData.bucket(node);
  const bool entityIsOwned = b.owned();
  const bool entityIsShared = b.shared();
  const bool entityIsGhosted = !entityIsOwned && !entityIsShared;

  if (realm_.hasPeriodic_ && realm_.has_non_matching_boundary_face_alg())
    throw std::logic_error("not ready for primetime to compbine periodic and non-matching algorithm suite");

  // simple case
  if (!realm_.hasPeriodic_ && !realm_.has_non_matching_boundary_face_alg()) {
    if (entityIsGhosted)
      return DS_GhostedDOF;
    if (entityIsOwned)
      return DS_OwnedDOF;
    if (entityIsShared && !entityIsOwned)
      return DS_GloballyOwnedDOF;
  }

  if (realm_.has_non_matching_boundary_face_alg()) {
    if (entityIsOwned)
      return DS_OwnedDOF;
    if (!entityIsOwned && (entityIsGhosted || entityIsShared)){
      return DS_GloballyOwnedDOF;
    }
  }

  if (realm_.hasPeriodic_) {
    const stk::mesh::EntityId stkId = bulkData.identifier(node);
    const stk::mesh::EntityId naluId = *stk::mesh::field_data(*realm_.naluGlobalId_, node);

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
  const unsigned p_size = bulkData.parallel_size();
  (void)p_size;

  // create a localID for all active nodes in the mesh...
  const stk::mesh::Selector s_universal = metaData.universal_part()
      & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& buckets =
      realm_.get_buckets( stk::topology::NODE_RANK, s_universal );

  // we allow for ghosted nodes when contact is active. When periodic is active, we may
  // also have ghosted nodes due to the periodicGhosting. However, we want to exclude these
  // nodes

  LocalOrdinal numGhostNodes = 0;
  LocalOrdinal numOwnedNodes = 0;
  LocalOrdinal numNodes = 0;
  LocalOrdinal numGloballyOwnedNotLocallyOwned = 0; // these are nodes on other procs
  // First, get the number of owned and globallyOwned (or num_globallyOwned_nodes = num_nodes - num_owned_nodes)
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ;
        ib != buckets.end() ; ++ib ) {
    const stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length = b.size();

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
  }

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
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ; ib != buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    {
      const stk::mesh::Bucket::size_type length   = b.size();
      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        const stk::mesh::Entity entity = b[k];

        int status = getDofStatus(entity);
        if (!(status & DS_SkippedDOF) && (status & DS_OwnedDOF))
          owned_nodes.push_back(entity);
      }
    }
  }
  
  std::sort(owned_nodes.begin(), owned_nodes.end(), CompareEntityById(bulkData, realm_.naluGlobalId_) );
  
  myLIDs_.clear();
  for (unsigned inode=0; inode < owned_nodes.size(); ++inode) {
    const stk::mesh::Entity entity = owned_nodes[inode];
    const stk::mesh::EntityId entityId = *stk::mesh::field_data(*realm_.naluGlobalId_, entity);
    // entityId can be duplicated in periodic or contact
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
  
  // now globallyOwned:
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ; ib != buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    {
      const stk::mesh::Bucket::size_type length   = b.size();
      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        const stk::mesh::Entity node = b[k];
        
        int status = getDofStatus(node);
        if (!(status & DS_SkippedDOF) && (status & DS_GloballyOwnedDOF))
          globally_owned_nodes.push_back(node);
      }
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

  const Teuchos::RCP<LinSys::Comm> tpetraComm = Tpetra::rcp(new LinSys::Comm(bulkData.parallel()));
  ownedRowsMap_ = Teuchos::rcp(new LinSys::Map(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), ownedGids, 1, tpetraComm, node_));
  globallyOwnedRowsMap_ = Teuchos::rcp(new LinSys::Map(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), globallyOwnedGids, 1, tpetraComm, node_));

  exporter_ = Teuchos::rcp(new LinSys::Export(globallyOwnedRowsMap_, ownedRowsMap_));
  importer_ = Teuchos::rcp(new LinSys::Import(ownedRowsMap_, globallyOwnedRowsMap_));

  ownedPlusGloballyOwnedRowsMap_ = Teuchos::rcp(new LinSys::Map(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), totalGids_, 1, tpetraComm, node_));

  globallyOwnedGraph_ = Teuchos::rcp(new LinSys::Graph(globallyOwnedRowsMap_, ownedPlusGloballyOwnedRowsMap_, 8));

  // Now, we're ready to have Algs call the build*Graph() methods and build up the connection list (row,col).
  // We'll finish this off in finalizeLinearSystem()
}

void TpetraLinearSystem::addConnections(const std::vector<stk::mesh::Entity> & entities)
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const unsigned p_rank = bulkData.parallel_rank();
  (void)p_rank;

  const size_t num_entities = entities.size();
  for(size_t a=0; a < num_entities; ++a) {
    const stk::mesh::Entity entity_a = entities[a];
    const stk::mesh::EntityId id_a = *stk::mesh::field_data(*realm_.naluGlobalId_, entity_a);

    for(size_t b=0; b < num_entities; ++b) {
      const stk::mesh::Entity entity_b = entities[b];
      const stk::mesh::EntityId id_b = *stk::mesh::field_data(*realm_.naluGlobalId_, entity_b);
      const bool a_then_b = id_a < id_b;
      const stk::mesh::Entity entity_min = a_then_b ? entity_a : entity_b;
      const stk::mesh::Entity entity_max = a_then_b ? entity_b : entity_a;
      connectionSet_.insert( Connection(entity_min, entity_max) );
    }
  }
}

void
TpetraLinearSystem::buildNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & metaData = realm_.meta_data();

  const stk::mesh::Selector s_owned = metaData.locally_owned_part()
    & stk::mesh::selectUnion(parts) 
    & !(stk::mesh::selectUnion(realm_.get_slave_part_vector()))
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_owned );
  std::vector<stk::mesh::Entity> entities(1);
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ;
        ib != buckets.end() ; ++ib ) {
    const stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      entities[0] = b[k];
      addConnections(entities);
    }
  }
}

void
TpetraLinearSystem::buildEdgeToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & metaData = realm_.meta_data();

  const stk::mesh::Selector s_owned = metaData.locally_owned_part()
    & stk::mesh::selectUnion(parts) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::EDGE_RANK, s_owned );
  const size_t numNodes = 2; // Edges are easy...
  std::vector<stk::mesh::Entity> entities(numNodes);
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ;
        ib != buckets.end() ; ++ib ) {
    const stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      stk::mesh::Entity const * edge_nodes = b.begin_nodes(k);

      // figure out the global dof ids for each dof on each node
      for(size_t n=0; n < numNodes; ++n) {
        entities[n] = edge_nodes[n];
      }
      addConnections(entities);
    }
  }
}

void
TpetraLinearSystem::buildFaceToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & metaData = realm_.meta_data();

  const stk::mesh::Selector s_owned = metaData.locally_owned_part()
    & stk::mesh::selectUnion(parts) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( metaData.side_rank(), s_owned );
  std::vector<stk::mesh::Entity> entities;
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ;
        ib != buckets.end() ; ++ib ) {
    const stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      stk::mesh::Entity const * face_nodes = b.begin_nodes(k);

      // figure out the global dof ids for each dof on each node
      const size_t numNodes = b.num_nodes(k);
      entities.resize(numNodes);
      for(size_t n=0; n < numNodes; ++n) {
        entities[n] = face_nodes[n];
      }
      addConnections(entities);
    }
  }
}

void
TpetraLinearSystem::buildElemToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & metaData = realm_.meta_data();

  const stk::mesh::Selector s_owned = metaData.locally_owned_part()
    & stk::mesh::selectUnion(parts) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_owned );
  std::vector<stk::mesh::Entity> entities;
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ;
        ib != buckets.end(); ++ib ) {
    const stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      stk::mesh::Entity const * elem_nodes = b.begin_nodes(k);
      // figure out the global dof ids for each dof on each node
      const size_t numNodes = b.num_nodes(k);
      entities.resize(numNodes);
      for(size_t n=0; n < numNodes; ++n) {
        entities[n] = elem_nodes[n];
      }
      addConnections(entities);
    }
  }
}

void
TpetraLinearSystem::buildReducedElemToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::MetaData & metaData = realm_.meta_data();

  const stk::mesh::Selector s_owned = metaData.locally_owned_part()
    & stk::mesh::selectUnion(parts) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_owned );
  std::vector<stk::mesh::Entity> entities;
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin() ;
        ib != buckets.end() ; ++ib ) {
    const stk::mesh::Bucket & b = **ib ;

    // extract master element
    MasterElement *meSCS = realm_.get_surface_master_element(b.topology());
    // extract master element specifics
    const int numScsIp = meSCS->numIntPoints_;
    const int *lrscv = meSCS->adjacentNodes();

    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      stk::mesh::Entity const * elem_nodes = b.begin_nodes(k);

      // figure out the global dof ids for each dof on each node
      const size_t numNodes = 2;
      entities.resize(numNodes);
      for (int j = 0; j < numScsIp; ++j){
        for(size_t n=0; n < numNodes; ++n) {
          entities[n] = elem_nodes[lrscv[2*j+n]];
        }
        addConnections(entities);
      }
    }
  }
}

void
TpetraLinearSystem::buildFaceElemToNodeGraph(const stk::mesh::PartVector & parts)
{
  beginLinearSystemConstruction();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  stk::mesh::MetaData & metaData = realm_.meta_data();

  const stk::mesh::Selector s_owned = metaData.locally_owned_part()
    & stk::mesh::selectUnion(parts) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( metaData.side_rank(), s_owned );
  std::vector<stk::mesh::Entity> entities;
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin() ;
        ib != face_buckets.end() ; ++ib ) {
    const stk::mesh::Bucket & b = **ib ;
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
      for(size_t n=0; n < numNodes; ++n) {
        entities[n] = elem_nodes[n];
      }
      addConnections(entities);
    }
  }
}

void
TpetraLinearSystem::buildEdgeHaloNodeGraph(
  const stk::mesh::PartVector &/*parts*/)
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  beginLinearSystemConstruction();

  std::vector<stk::mesh::Entity> entities;

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
      stk::mesh::Entity const* elem_nodes = bulkData.begin_nodes(elem);
      const size_t numNodes = bulkData.num_nodes(elem);
      const size_t numEntities = numNodes+1;
      entities.resize(numEntities);

      entities[0] = node;
      for(size_t n=0; n < numNodes; ++n) {
        entities[n+1] = elem_nodes[n];
      }
      addConnections(entities);
    }
  }
}

void
TpetraLinearSystem::buildNonConformalNodeGraph(
  const stk::mesh::PartVector &/*parts*/)
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  beginLinearSystemConstruction();

  std::vector<stk::mesh::Entity> entities;

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
        for ( int ni = 0; ni < current_num_elem_nodes; ++ni ) {
          entities[ni] = current_elem_node_rels[ni];
        }
        
        // fill in connected nodes; opposing
        for ( int ni = 0; ni < opposing_num_elem_nodes; ++ni ) {
          entities[current_num_elem_nodes+ni] = opposing_elem_node_rels[ni];
        }
        
        // okay, now add the connections; will be symmetric 
        // columns of current node row (opposing nodes) will add columns to opposing nodes row
        addConnections(entities);
      }
    }
  }
}

void
TpetraLinearSystem::buildOversetNodeGraph(
  const stk::mesh::PartVector &/*parts*/)
{
  // extract the rank
  const int theRank = NaluEnv::self().parallel_rank();

  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  beginLinearSystemConstruction();

  std::vector<stk::mesh::Entity> entities;

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
    addConnections(entities);
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

  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin();
        ib != buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

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
TpetraLinearSystem::finalizeLinearSystem()
{
  ThrowRequire(inConstruction_);
  inConstruction_ = false;

  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  stk::mesh::MetaData & metaData = realm_.meta_data();

  const int this_mpi_rank = bulkData.parallel_rank();
  (void)this_mpi_rank;

  ConnectionVec connectionVec(connectionSet_.begin(), connectionSet_.end());
  connectionSet_.clear();
  std::sort(connectionVec.begin(), connectionVec.end());

  std::vector<GlobalOrdinal> globalDofs_a(numDof_);
  std::vector<GlobalOrdinal> globalDofs_b(numDof_);
  std::ostringstream out2;
  const size_t numConnections = connectionVec.size();
  for (size_t i=0; i < numConnections; ++i) {
    const stk::mesh::Entity entity_a = connectionVec[i].first;
    const stk::mesh::Entity entity_b = connectionVec[i].second;

    const stk::mesh::EntityId entityId_a = *stk::mesh::field_data(*realm_.naluGlobalId_, entity_a);
    const stk::mesh::EntityId entityId_b = *stk::mesh::field_data(*realm_.naluGlobalId_, entity_b);

    for (size_t d=0; d < numDof_; ++d) {
      globalDofs_a[d] = GID_(entityId_a, numDof_ , d);
      globalDofs_b[d] = GID_(entityId_b, numDof_ , d);
    }

    // NOTE: 'Connections' should already include the self
    // pairings (where entity_a == entity_b) so we don't have
    // to worry about doing an insert on (globalRow_a, globalDofs_a),
    // etc.

    // for dofs on entity_a add columns due to entity_b dofs
    if (getDofStatus(entity_a) & DS_GloballyOwnedDOF) { // !Locally owned
      for (size_t d=0; d < numDof_; ++d) {
        const GlobalOrdinal globalRow_a = GID_(entityId_a, numDof_, d);
        globallyOwnedGraph_->insertGlobalIndices(globalRow_a, globalDofs_b);
      }
    }

    // for dofs on entity_b add columns due to entity_a dofs
    if (getDofStatus(entity_b) & DS_GloballyOwnedDOF) { // !Locally owned
      for (size_t d=0; d < numDof_; ++d) {
        const GlobalOrdinal globalRow_b = GID_(entityId_b, numDof_ , d);
        globallyOwnedGraph_->insertGlobalIndices(globalRow_b, globalDofs_a);
      }
    }
  }

  std::ostringstream out3;
  globallyOwnedGraph_->fillComplete();

  LinSys::Graph ownedPlusGloballyOwnedGraph(ownedRowsMap_, 8);
  ownedPlusGloballyOwnedGraph.doExport(*globallyOwnedGraph_, *exporter_, Tpetra::INSERT);
  ownedPlusGloballyOwnedGraph.fillComplete(ownedRowsMap_, ownedRowsMap_);

  // Add columns that are imported to the totalGids_ array
  const Teuchos::RCP<const LinSys::Map> & map = ownedPlusGloballyOwnedGraph.getColMap();
  for (size_t i=0; i < map->getNodeNumElements(); ++i) {
    const GlobalOrdinal gid = map->getGlobalElement(i);
    if(!ownedPlusGloballyOwnedRowsMap_->isNodeGlobalElement(gid))
      totalGids_.push_back(gid);
  }

  // This is the column map for the owned graph now
  const Teuchos::RCP<LinSys::Comm> tpetraComm = Tpetra::rcp(new LinSys::Comm(bulkData.parallel()));
  totalColsMap_ = Teuchos::rcp(new LinSys::Map(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), totalGids_, 1, tpetraComm, node_));
  ownedGraph_ = Teuchos::rcp(new LinSys::Graph(ownedRowsMap_, totalColsMap_, 8));

  // Insert all the local connection data
  for (size_t i=0; i < numConnections; ++i) {
    const stk::mesh::Entity entity_a = connectionVec[i].first;
    const stk::mesh::Entity entity_b = connectionVec[i].second;

    const stk::mesh::EntityId entityId_a = *stk::mesh::field_data(*realm_.naluGlobalId_, entity_a);
    const stk::mesh::EntityId entityId_b = *stk::mesh::field_data(*realm_.naluGlobalId_, entity_b);

    for (size_t d=0; d < numDof_; ++d) {
      globalDofs_a[d] = GID_(entityId_a, numDof_, d);
      globalDofs_b[d] = GID_(entityId_b, numDof_, d);
    }

    // NOTE: 'Connections' should already include the self
    // pairings (where entity_a == entity_b) so we don't have
    // to worry about doing an insert on (globalRow_a, globalDofs_a),
    // etc.

    if (getDofStatus(entity_a) & DS_OwnedDOF) { // Locally owned
      for (size_t d=0; d < numDof_; ++d) {
        const GlobalOrdinal globalRow_a = GID_(entityId_a, numDof_ , d);
        ownedGraph_->insertGlobalIndices(globalRow_a, globalDofs_b);
      }
    }

    if (getDofStatus(entity_b) & DS_OwnedDOF) { // Locally owned
      for (size_t d=0; d < numDof_; ++d) {
        const GlobalOrdinal globalRow_b = GID_(entityId_b, numDof_ , d);
        ownedGraph_->insertGlobalIndices(globalRow_b, globalDofs_a);
      }
    }
  }

  // add imported graph information
  {
    const LinSys::Map & rowMap = *ownedPlusGloballyOwnedGraph.getRowMap();
    const LinSys::Map & colMap = *ownedPlusGloballyOwnedGraph.getColMap();
    const size_t numRows = rowMap.getNodeNumElements();
    std::vector<GlobalOrdinal> newInd;
    for(size_t localRow=0; localRow < numRows; ++localRow) {
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

  ownedRhs_ = Teuchos::rcp(new LinSys::Vector(ownedRowsMap_));
  globallyOwnedRhs_ = Teuchos::rcp(new LinSys::Vector(globallyOwnedRowsMap_));

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
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  const size_t n_obj = entities.size();
  const size_t numRows = n_obj * numDof_;

  ThrowAssert(numRows == rhs.size());
  ThrowAssert(numRows*numRows == lhs.size());

  for(size_t i=0; i < n_obj; ++i) {
    const stk::mesh::Entity entity = entities[i];
    const stk::mesh::EntityId entityId = bulkData.identifier(entity);
    (void)entityId;
    const stk::mesh::EntityId naluId = *stk::mesh::field_data(*realm_.naluGlobalId_, entity);
    const LocalOrdinal localOffset = lookup_myLID(myLIDs_, naluId, "sumInto", entity) * numDof_;
    for(size_t d=0; d < numDof_; ++d) {
      size_t lid = i*numDof_ + d;
      scratchIds[lid] = localOffset + d;
    }
  }

  for(size_t r=0; r < numRows; ++r) {
    const LocalOrdinal localId = scratchIds[r];

    for(size_t c=0; c < numRows; ++c) // numRows == numCols
      scratchVals[c] = lhs[r*numRows + c];

    if(localId < maxOwnedRowId_) {
      ownedMatrix_->sumIntoLocalValues(localId, scratchIds, scratchVals);
      ownedRhs_->sumIntoLocalValue(localId, rhs[r]);
    }
    else if(localId < maxGloballyOwnedRowId_) {
      const LocalOrdinal actualLocalId = localId - maxOwnedRowId_;
      globallyOwnedMatrix_->sumIntoLocalValues(actualLocalId, scratchIds, scratchVals);
      globallyOwnedRhs_->sumIntoLocalValue(actualLocalId, rhs[r]);
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
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  double adbc_time = -stk::cpu_time();
  const unsigned p_size = bulkData.parallel_size();
  (void)p_size;

  const stk::mesh::Selector selector 
    = (metaData.locally_owned_part() | metaData.globally_shared_part())
    & stk::mesh::selectUnion(parts)
    & stk::mesh::selectField(*solutionField) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, selector );

  int nbc=0;
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin();
        ib != buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

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

        if(localId > maxGloballyOwnedRowId_) {
          std::cout << "localId > maxGloballyOwnedRowId_:: localId= " << localId << " maxGloballyOwnedRowId_= " << maxGloballyOwnedRowId_ << " useOwned = " << (localId < maxOwnedRowId_ ) << std::endl;
          throw std::runtime_error("logic error: localId > maxGloballyOwnedRowId_");
        }

        // Adjust the LHS

        const double diagonal_value = useOwned ? 1.0 : 0.0;

        matrix->getLocalRowView(actualLocalId, indices, values);
        const size_t rowLength = values.size();
        new_values.resize(rowLength);
        for(size_t i=0; i < rowLength; ++i) {
            new_values[i] = (indices[i] == localId) ? diagonal_value : 0;
        }
        matrix->replaceLocalValues(actualLocalId, indices, new_values);

        // Replace the RHS residual with (desired - actual)
        Teuchos::RCP<LinSys::Vector> rhs = useOwned ? ownedRhs_: globallyOwnedRhs_;
        const double bc_residual = useOwned ? (bcValues[k*fieldSize + d] - solution[k*fieldSize + d]) : 0.0;
        rhs->replaceLocalValue(actualLocalId, bc_residual);
        ++nbc;
      }
    }
  }
  adbc_time += stk::cpu_time();
}

void
TpetraLinearSystem::prepareConstraints(
  const unsigned beginPos,
  const unsigned endPos)
{
  Teuchos::ArrayView<const LocalOrdinal> indices;
  Teuchos::ArrayView<const double> values;
  std::vector<double> new_values;

  // iterate oversetInfoVec_
  std::vector<OversetInfo *>::iterator ii;
  for( ii=realm_.oversetManager_->oversetInfoVec_.begin();
       ii!=realm_.oversetManager_->oversetInfoVec_.end(); ++ii ) {

    // extract orphan node and global id; process both owned and shared
    stk::mesh::Entity orphanNode = (*ii)->orphanNode_;
    const stk::mesh::EntityId naluId = *stk::mesh::field_data(*realm_.naluGlobalId_, orphanNode);
    const LocalOrdinal localIdOffset = lookup_myLID(myLIDs_, naluId, "prepareConstraints") * numDof_;

    for(unsigned d=beginPos; d < endPos; ++d) {
      const LocalOrdinal localId = localIdOffset + d;
      const bool useOwned = localId < maxOwnedRowId_;
      const LocalOrdinal actualLocalId = useOwned ? localId : localId - maxOwnedRowId_;
      Teuchos::RCP<LinSys::Matrix> matrix = useOwned ? ownedMatrix_ : globallyOwnedMatrix_;
      
      if ( localId > maxGloballyOwnedRowId_) {
        throw std::runtime_error("logic error: localId > maxGloballyOwnedRowId_");
      }
      
      // Adjust the LHS; full row is perfectly zero
      matrix->getLocalRowView(actualLocalId, indices, values);
      const size_t rowLength = values.size();
      new_values.resize(rowLength);
      for(size_t i=0; i < rowLength; ++i) {
        new_values[i] = 0.0;
      }
      matrix->replaceLocalValues(actualLocalId, indices, new_values);
      
      // Replace the RHS residual with zero
      Teuchos::RCP<LinSys::Vector> rhs = useOwned ? ownedRhs_: globallyOwnedRhs_;
      const double bc_residual = 0.0;
      rhs->replaceLocalValue(actualLocalId, bc_residual);
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
    writeToFile(this->name_.c_str());
    writeToFile(this->name_.c_str(), false);
  }

  double solve_time = -stk::cpu_time();

  int iters;
  double finalResidNorm;
  
  // memory diagnostic
  if ( realm_.get_activate_memory_diagnostic() ) {
    NaluEnv::self().naluOutputP0() << "NaluMemory::TpetraLinearSystem::solve() PreSolve: " << name_ << std::endl;
    realm_.provide_memory_summary();
  }

  const int status = linearSolver->solve(
      sln_,
      iters,
      finalResidNorm);

  solve_time += stk::cpu_time();

  if (linearSolver->getConfig()->getWriteMatrixFiles()) {
    writeSolutionToFile(this->name_.c_str());
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
TpetraLinearSystem::checkForNaN(bool useOwned)
{
  Teuchos::RCP<LinSys::Matrix> matrix = useOwned ? ownedMatrix_ : globallyOwnedMatrix_;
  Teuchos::RCP<LinSys::Vector> rhs = useOwned ? ownedRhs_ : globallyOwnedRhs_;

  Teuchos::ArrayView<const LocalOrdinal> indices;
  Teuchos::ArrayView<const double> values;

  int n = matrix->getRowMap()->getNodeNumElements();
  for (int i=0; i < n; ++i) {
    matrix->getLocalRowView(i, indices, values);
    const size_t rowLength = values.size();
    for(size_t k=0; k < rowLength; ++k) {
      if (values[k] != values[k])	{
        std::cout << "LHS NaN: " << i << std::endl;
        throw std::runtime_error("bad LHS");
      }
    }
  }

  Teuchos::ArrayRCP<const Scalar> rhs_data = rhs->getData();
  n = rhs_data.size();
  for (int i=0; i < n; ++i) {
    if (rhs_data[i] != rhs_data[i]) {
      std::cout << "rhs NaN: " << i << std::endl;
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
  int n = matrix->getRowMap()->getNodeNumElements();
  GlobalOrdinal max_gid = 0, g_max_gid=0;
  for (int i=0; i < n; ++i) {
    GlobalOrdinal gid = matrix->getGraph()->getRowMap()->getGlobalElement(i);
    max_gid = std::max(gid, max_gid);
  }
  stk::all_reduce_max(bulkData.parallel(), &max_gid, &g_max_gid, 1);

  nrowG = g_max_gid+1;
  std::vector<double> local_row_sums(nrowG, 0.0);
  std::vector<int> local_row_exists(nrowG, 0);
  std::vector<double> global_row_sums(nrowG, 0.0);
  std::vector<int> global_row_exists(nrowG, 0);

  for (int i=0; i < n; ++i) {
    GlobalOrdinal gid = matrix->getGraph()->getRowMap()->getGlobalElement(i);
    matrix->getLocalRowView(i, indices, values);
    const size_t rowLength = values.size();
    double row_sum = 0.0;
    for(size_t k=0; k < rowLength; ++k) {
      row_sum += std::abs(values[k]);
    }
    if (gid-1 >= (GlobalOrdinal)local_row_sums.size() || gid <= 0) {
      std::cout << "gid= " << gid << " nrowG= " << nrowG << std::endl;
      throw std::runtime_error("bad gid");
    }
    local_row_sums[gid-1] = row_sum;
    local_row_exists[gid-1] = 1;
  }

  stk::all_reduce_sum(bulkData.parallel(), &local_row_sums[0], &global_row_sums[0], (unsigned)nrowG);
  stk::all_reduce_max(bulkData.parallel(), &local_row_exists[0], &global_row_exists[0], (unsigned)nrowG);

  bool found=false;
  for (size_t ii=0; ii < nrowG; ++ii) {
    double row_sum = global_row_sums[ii];
    if (global_row_exists[ii] && bulkData.parallel_rank() == 0 && row_sum < 1.e-10) {
      found = true;
      GlobalOrdinal gid = ii+1;
      stk::mesh::EntityId nid = GLOBAL_ENTITY_ID(gid, numDof_);
      std::cout << "nid= " << nid << std::endl;
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
  }

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
                                                                    name_, std::string("Tpetra matrix for: ")+name_, true);
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

  if (p_rank == 0)
    {
      std::cout << "\nMatrix for system: " << name_ << " :: N N NZ= " << matrix->getRangeMap()->getGlobalNumElements()
                << " "
                << matrix->getDomainMap()->getGlobalNumElements()
                << " "
                << matrix->getGlobalNumEntries()
                << std::endl;
      NaluEnv::self().naluOutputP0() << "\nMatrix for system: " << name_ << " :: N N NZ= " << matrix->getRangeMap()->getGlobalNumElements()
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
  (void)p_rank;

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
      const LocalOrdinal localIdOffset = lookup_myLID(myLIDs_, naluGlobalId[k], "copy_tpetra_to_stk") * numDof_;
      stk::mesh::Entity node = b[k];
      stk::mesh::EntityId stkId = bulkData.identifier(node);
      stk::mesh::EntityId naluId = naluGlobalId[k];
      for(unsigned d=0; d < fieldSize; ++d) {
        const LocalOrdinal localId = localIdOffset + d;
        bool useOwned = true;
        LocalOrdinal actualLocalId = localId;
        if(localId >= maxOwnedRowId_) {
          actualLocalId = localId - maxOwnedRowId_;
          useOwned = false;
        }

        if (!useOwned) {
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

} // namespace nalu
} // namespace Sierra
