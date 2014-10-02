/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <PeriodicManager.h>
#include <NaluEnv.h>
#include <Realm.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_io
#include <stk_io/StkMeshIoBroker.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/CPUTime.hpp>

// stk_search
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>

// vector
#include <vector>
#include <map>

namespace sierra{
namespace nalu{

PeriodicManager::PeriodicManager(
   Realm &realm)
  : realm_(realm ),
    searchTolerance_(1.0e8),
    periodicGhosting_(NULL),
    ghostingName_("nalu_periodic"),
    timerSearch_(0.0)
{
  // do nothing
}

PeriodicManager::~PeriodicManager()
{
  // nothing to delete
}

//--------------------------------------------------------------------------
//-------- add periodic pair -----------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::add_periodic_pair(
  stk::mesh::Part * masterMeshPart,
  stk::mesh::Part * slaveMeshPart,
  const double &searchTolerance,
  const std::string &searchMethodName)
{
  // use most stringent tolerance
  if (searchTolerance < searchTolerance_)
    searchTolerance_ = searchTolerance;

  MeshPartPair periodicPartPair(masterMeshPart, slaveMeshPart);
  periodicPartPairs_.push_back(periodicPartPair);

  // determine search method for this pair; default is boost_rtree
  stk::search::SearchMethod searchMethod = stk::search::BOOST_RTREE;
  if ( searchMethodName == "boost_rtree" )
      searchMethod = stk::search::BOOST_RTREE;
  else if ( searchMethodName == "stk_octree" )
    searchMethod = stk::search::OCTREE;
  else
    NaluEnv::self().naluOutputP0() << "PeriodicManager::search method not declared; will use BOOST_RTREE" << std::endl;
  searchMethodVec_.push_back(searchMethod);
}

//--------------------------------------------------------------------------
//-------- get_slave_part_vector ----------------------------------------------
//--------------------------------------------------------------------------
const stk::mesh::PartVector &
PeriodicManager::get_slave_part_vector()
{
  return slavePartVector_;
}

//--------------------------------------------------------------------------
//-------- get_search_time -------------------------------------------------
//--------------------------------------------------------------------------
double
PeriodicManager::get_search_time() {
  return timerSearch_;
}

//--------------------------------------------------------------------------
//-------- build_constraints -----------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::build_constraints()
{

  if ( periodicPartPairs_.size() == 0 )
    throw std::runtime_error("PeriodiocBC::Error: No periodic pair provided");

  // multiple pairs supported, however, constraint resolution not..
  if ( periodicPartPairs_.size() > 1 ) {
    NaluEnv::self().naluOutputP0()
      << "Multiple periodic pairs supported, however, must not have any reduction need" << std::endl;
  }

  // translate, search, constraint mapping
  for ( size_t k = 0; k < periodicPartPairs_.size(); ++k) {
    setup_gid_pairs(periodicPartPairs_[k].first, periodicPartPairs_[k].second, searchMethodVec_[k]);
  }

  create_ghosting_object();

  master_slave_reduction();

  create_slave_part_vector();

  // provide Nalu id update
  update_global_id_field();

}

//--------------------------------------------------------------------------
//-------- setup_gid_pairs -------------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::setup_gid_pairs(
  const stk::mesh::Part *masterPart,
  const stk::mesh::Part *slavePart,
  const stk::search::SearchMethod searchMethod)
{

  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "Periodic Review:  realm: " << realm_.name_ << std::endl;
  NaluEnv::self().naluOutputP0() << "=========================" << std::endl;

  const std::string masterPartName = masterPart->name();
  const std::string slavePartName = slavePart->name();

  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();
  stk::mesh::BulkData & bulk_data = realm_.fixture_->bulk_data();

  // required data structures; master/slave
  std::vector<sphereBoundingBox> sphereBoundingBoxMasterVec;
  std::vector<sphereBoundingBox> sphereBoundingBoxSlaveVec;

  // fields
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  const int nDim = meta_data.spatial_dimension();

  // hack due to 3d; 2d searches need to fake a 3d search
  const double pointRadius = searchTolerance_;
  Point masterCenter, slaveCenter;

  // Master: global_sum_coords_master; setup sphereBoundingBoxMasterVec,
  std::vector<double> local_sum_coords_master(nDim, 0.0), global_sum_coords_master(nDim, 0.0);

  size_t numberMasterNodes = 0, g_numberMasterNodes = 0;

  stk::mesh::Selector sm_locally_owned = meta_data.locally_owned_part()
    &stk::mesh::Selector(*masterPart);
  stk::mesh::BucketVector const& master_node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, sm_locally_owned );

  for ( stk::mesh::BucketVector::const_iterator ib = master_node_buckets.begin();
        ib != master_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length   = b.size();
    // point to data
    const double * coords = stk::mesh::field_data(*coordinates, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      stk::mesh::Entity node = b[k];
      numberMasterNodes += 1;
      // setup ident
      theEntityKey theIdent(bulk_data.entity_key(node), NaluEnv::self().parallel_rank());

      // define offset for all nodal fields that are of nDim
      const size_t offSet = k*nDim;

      // sum local coords for translation; define localCoords for bounding point
      for (int j = 0; j < nDim; ++j ) {
        const double cxj = coords[offSet+j];
        local_sum_coords_master[j] += cxj;
        masterCenter[j] = cxj;
      }

      // create the bounding point sphere and push back
      sphereBoundingBox theSphere( Sphere(masterCenter, pointRadius), theIdent);
      sphereBoundingBoxMasterVec.push_back(theSphere);
    }
  }
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &local_sum_coords_master[0], &global_sum_coords_master[0], nDim);
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &numberMasterNodes, &g_numberMasterNodes, 1);

  // Slave: global_sum_coords_slave
  std::vector<double> local_sum_coords_slave(nDim, 0.0), global_sum_coords_slave(nDim, 0.0);

  stk::mesh::Selector ss_locally_owned = meta_data.locally_owned_part()
    &stk::mesh::Selector(*slavePart);
  stk::mesh::BucketVector const& slave_node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, ss_locally_owned );

  size_t numberSlaveNodes = 0, g_numberSlaveNodes = 0;

  for ( stk::mesh::BucketVector::const_iterator ib = slave_node_buckets.begin();
        ib != slave_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length   = b.size();
    // point to data
    const double * coords = stk::mesh::field_data(*coordinates, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      numberSlaveNodes += 1;
      // define offset for all nodal fields that are of nDim
      const size_t offSet = k*nDim;
      for (int j = 0; j < nDim; ++j ) {
        local_sum_coords_slave[j] += coords[offSet+j];
      }
    }
  }
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &local_sum_coords_slave[0], &global_sum_coords_slave[0], nDim);
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &numberSlaveNodes, &g_numberSlaveNodes, 1);

  // throw right away if master and slave nodes are not equal
  if ( g_numberMasterNodes != g_numberSlaveNodes ) {
    NaluEnv::self().naluOutputP0() << "Periodic BC Error: Master/Slave Pair: "
        << masterPartName << "/" << slavePartName << std::endl;
    NaluEnv::self().naluOutputP0() << "Master part has " << g_numberMasterNodes
        << " while Slave part has " << g_numberSlaveNodes << std::endl;
    throw std::runtime_error("Please ensure that periodic part surfaces match - both in node count and sane connectivities");
  }

  // TRANSLATE SLAVE onto MASTER
  std::vector<double> translationVector(nDim, 0.0);
  std::vector<double> rotationPoint(nDim, 0.0);  // master center of mass

  for (int j = 0; j < nDim; ++j ) {
    translationVector[j] = (global_sum_coords_master[j] - global_sum_coords_slave[j]) / g_numberMasterNodes;
    rotationPoint[j] = global_sum_coords_master[j] / g_numberMasterNodes;
  }

  NaluEnv::self().naluOutputP0() << "Translating [ ";
  for (int j = 0; j < nDim; ++j ) {  NaluEnv::self().naluOutputP0() << translationVector[j] << " "; }
  NaluEnv::self().naluOutputP0() << "] for Master/Slave pair " << masterPartName << "/" << slavePartName << std::endl;

  // SLAVE - setup sphereBoundingBoxSlaveVec
  for ( stk::mesh::BucketVector::const_iterator ib = slave_node_buckets.begin();
       ib != slave_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length   = b.size();
    // point to data
    const double * coords = stk::mesh::field_data(*coordinates, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      stk::mesh::Entity node = b[k];
      // setup ident
      theEntityKey theIdent(bulk_data.entity_key(node), NaluEnv::self().parallel_rank());
      // define offset for all nodal fields that are of nDim
      const size_t offSet = k*nDim;

      // create the bounding point that is translated, then push back
      for (int j = 0; j < nDim; ++j ) {
        const double xj = coords[offSet+j];
        slaveCenter[j] = xj + translationVector[j];
      }
      sphereBoundingBox theSphere(Sphere(slaveCenter, pointRadius), theIdent);
      sphereBoundingBoxSlaveVec.push_back(theSphere);
    }
  }

  // will want to stuff product of search to a single vector
  std::vector<std::pair<theEntityKey, theEntityKey> > searchKeyPair;
  double timeA = stk::cpu_time();
  stk::search::coarse_search(sphereBoundingBoxSlaveVec, sphereBoundingBoxMasterVec, searchMethod, NaluEnv::self().parallel_comm(), searchKeyPair);
  timerSearch_ += (stk::cpu_time() - timeA);

  //=====================================================================
  // begin a series of checks that might indicate an issue with the mesh
  //=====================================================================

  // search should have provided **something**
  if ( searchKeyPair.size() < 1 )
    Env::output() << " issue with coarse_search; no pairs found. try increasing tolerance " << std::endl;

  // Each slave entry must have one master. If not, smaller tolerance is required
  const size_t searchSize = searchKeyPair.size();
  size_t problemNodes = 0;
  std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId>  > problemPairVec;
  for ( size_t k = 0; k < searchSize; ++k) {
    // get slave node for kth entry
    stk::mesh::EntityId kthSlaveId = searchKeyPair[k].first.id();
    size_t maxSizeToCheck = searchSize - 1;
    if ( k < maxSizeToCheck ) {
      // okay to check next entry
      stk::mesh::EntityId kp1SlaveId = searchKeyPair[k+1].first.id();

      // are they the same? If so, we have a problem
      if ( kp1SlaveId == kthSlaveId ) {
        stk::mesh::EntityId kthMasterId = searchKeyPair[k].second.id();
        std::pair<stk::mesh::EntityId, stk::mesh::EntityId> theProblemPair
          = std::make_pair(kthMasterId, kthSlaveId);
        problemPairVec.push_back(theProblemPair);
        problemNodes++;
      }
    }
  }
  size_t g_problemNodes = 0;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &problemNodes, &g_problemNodes, 1);

  // report issues
  if ( g_problemNodes > 0 ){
    NaluEnv::self().naluOutputP0() << "PeriodicSearchError: Multiple candidates for periodic pair search on pair "
        << masterPartName <<"/" << slavePartName << std::endl << std::endl;
    NaluEnv::self().naluOutputP0() << g_problemNodes << " Occurrences; "
        << "Please reduce search_tolerance to " << searchTolerance_/10.0 << std::endl;
    for ( size_t k = 0; k < problemPairVec.size(); ++k ) {
      stk::mesh::EntityId masterId = problemPairVec[k].first;
      stk::mesh::EntityId slaveId = problemPairVec[k].second;
      NaluEnv::self().naluOutputP0() << "Candidate Master/Slave pairs " << masterId << "/" << slaveId << std::endl;
    }
    throw std::runtime_error("PeriodiocBC::Error: Please reduce periodic_search_tolerance");
  }

  // All of the master nodes must have found at least one slave node
  problemNodes = 0;
  std::vector<stk::mesh::EntityId> problemNodeVec;
  for ( stk::mesh::BucketVector::const_iterator ib = master_node_buckets.begin();
          ib != master_node_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib;
      const stk::mesh::Bucket::size_type length   = b.size();
      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

        stk::mesh::EntityId masterIdFromPart = bulk_data.identifier(b[k]);

        // loop over search and see if one exists
        bool foundIt = false;
        for ( size_t j = 0; j < searchSize; ++j) {
          stk::mesh::EntityId masterIdFromSearch = searchKeyPair[j].second.id();
          if ( masterIdFromPart == masterIdFromSearch ) {
            foundIt = true;
            break;
          }
        }
        // if we did not find the node, we have a problem
        if (!foundIt ) {
          problemNodes++;
          problemNodeVec.push_back(masterIdFromPart);
        }
      }
  }
  g_problemNodes = 0;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &problemNodes, &g_problemNodes, 1);

  // report issues
  if ( g_problemNodes > 0 ){
    NaluEnv::self().naluOutputP0() << "PeriodicSearchError: Master node never found a slave node match: "
        << masterPartName <<"/" << slavePartName << std::endl << std::endl;
    NaluEnv::self().naluOutputP0() << g_problemNodes << " Occurrences; "
         << "Please check your mesh for proper periodicity " << std::endl;
     for ( size_t k = 0; k < problemNodeVec.size(); ++k ) {
       NaluEnv::self().naluOutputP0() << "Problem Master node " << problemNodeVec[k] << std::endl;
     }
     throw std::runtime_error("PeriodiocBC::Error: Master node did not find a slave node; increase search tolerance?");
  }

  // populate searchKeyVector_; culmination of all master/slaves
  searchKeyVector_.insert(searchKeyVector_.end(), searchKeyPair.begin(), searchKeyPair.end());

}

//--------------------------------------------------------------------------
//-------- create_ghosting_object ------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::create_ghosting_object()
{

  stk::mesh::BulkData & bulk_data = realm_.fixture_->bulk_data();
  unsigned theRank = NaluEnv::self().parallel_rank();

  std::vector<stk::mesh::EntityProc> sendNodes;
  for (size_t i=0, size=searchKeyVector_.size(); i<size; ++i) {

    unsigned domainProc = searchKeyVector_[i].first.proc();
    unsigned rangeProc = searchKeyVector_[i].second.proc();

    if ((theRank != domainProc) && (theRank == rangeProc)) {
      stk::mesh::Entity rangeNode = bulk_data.get_entity(searchKeyVector_[i].second.id());
      sendNodes.push_back(stk::mesh::EntityProc(rangeNode, domainProc));
    }
    else if ((theRank == domainProc) && (theRank != rangeProc)) {
      stk::mesh::Entity domainNode = bulk_data.get_entity(searchKeyVector_[i].first.id());
      sendNodes.push_back(stk::mesh::EntityProc(domainNode, rangeProc));
    }
  }

  size_t numNodes = sendNodes.size();
  size_t g_numNodes = 0;
  stk::all_reduce_sum(bulk_data.parallel(), &numNodes, &g_numNodes, 1);
  if ( g_numNodes > 0) {
    // check if we need to ghost
    bulk_data.modification_begin();
    periodicGhosting_ = &bulk_data.create_ghosting(ghostingName_);
    bulk_data.change_ghosting(*periodicGhosting_, sendNodes);
    bulk_data.modification_end();
  }

  // now populate master slave communicator
  for (size_t i=0, size=searchKeyVector_.size(); i<size; ++i) {
     stk::mesh::Entity domainNode = bulk_data.get_entity(searchKeyVector_[i].first.id());
     stk::mesh::Entity rangeNode = bulk_data.get_entity(searchKeyVector_[i].second.id());
     // unique master:slave communicator
     EntityPair theFirstPair = std::make_pair(rangeNode, domainNode);
     masterSlaveCommunicator_.push_back(theFirstPair);
  }

}

//--------------------------------------------------------------------------
//-------- get_ghosting_object ------------------------------------------
//--------------------------------------------------------------------------
stk::mesh::Ghosting *
PeriodicManager::get_ghosting_object()
{
  return periodicGhosting_;
}

//--------------------------------------------------------------------------
//-------- master_slave_reduction ------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::master_slave_reduction()
{
  // multiple pairs supported, however, constraint resolution not..
  if ( periodicPartPairs_.size() > 1 ) {
    NaluEnv::self().naluOutputP0()
      << "Multiple periodic pairs supported, however, must not have any reduction need" << std::endl;
  }
}

//--------------------------------------------------------------------------
//-------- create_slave_part_vector ----------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::create_slave_part_vector()
{
  // for now, no slave reduction
  for ( size_t k = 0; k < periodicPartPairs_.size(); ++k) {
    slavePartVector_.push_back(periodicPartPairs_[k].second);
  }
}

//--------------------------------------------------------------------------
//-------- periodic_parallel_communicate_field -----------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::periodic_parallel_communicate_field(
  stk::mesh::FieldBase *theField)
{
  if ( NULL != periodicGhosting_ ) {
    std::vector< const stk::mesh::FieldBase *> fieldVec(1, theField);
    stk::mesh::communicate_field_data(*periodicGhosting_, fieldVec);
  }
}

//--------------------------------------------------------------------------
//-------- parallel_communicate_field --------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::parallel_communicate_field(
  stk::mesh::FieldBase *theField)
{
  stk::mesh::BulkData & bulk_data = realm_.fixture_->bulk_data();
  const unsigned pSize = bulk_data.parallel_size();
  if ( pSize > 1 ) {
    std::vector< const stk::mesh::FieldBase *> fieldVec(1, theField);
    stk::mesh::copy_owned_to_shared( bulk_data, fieldVec);
    stk::mesh::communicate_field_data(bulk_data.aura_ghosting(), fieldVec);
  }
}

//--------------------------------------------------------------------------
//-------- update global_id_field ------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::update_global_id_field()
{

  stk::mesh::BulkData & bulk_data = realm_.fixture_->bulk_data();

  // no need to update periodically ghosted fields

  // vector of masterEntity:slaveEntity pairs
  for ( size_t k = 0; k < masterSlaveCommunicator_.size(); ++k) {

    // extract master node and slave node
    EntityPair vecPair = masterSlaveCommunicator_[k];
    const stk::mesh::Entity masterNode = vecPair.first;
    const stk::mesh::Entity slaveNode = vecPair.second;

    // pointer to data
    const stk::mesh::EntityId masterGlobalId = bulk_data.identifier(masterNode );
    stk::mesh::EntityId *slaveGlobalId = stk::mesh::field_data(*realm_.naluGlobalId_, slaveNode);

    // set new value
    *slaveGlobalId = masterGlobalId;

  }

  // update all shared; aura and periodic
  parallel_communicate_field(realm_.naluGlobalId_);

}

//--------------------------------------------------------------------------
//-------- apply_constraints -----------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::apply_constraints(
  stk::mesh::FieldBase *theField,
  const unsigned &sizeOfField,
  const bool &bypassFieldCheck,
  const bool &addSlaves,
  const bool &setSlaves)
{

  // update periodically ghosted fields
  periodic_parallel_communicate_field(theField);

  if ( addSlaves )
    add_slave_to_master(theField, sizeOfField, bypassFieldCheck);
  if ( setSlaves )
    set_slave_to_master(theField, sizeOfField, bypassFieldCheck);

  // parallel communicate shared and aura-ed entities
  parallel_communicate_field(theField);

}


//--------------------------------------------------------------------------
//-------- apply_max_field -------------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::apply_max_field(
  stk::mesh::FieldBase *theField,
  const unsigned &sizeOfField)
{

  for ( size_t k = 0; k < masterSlaveCommunicator_.size(); ++k) {
    // extract master node and slave node
    EntityPair vecPair = masterSlaveCommunicator_[k];
    const stk::mesh::Entity masterNode = vecPair.first;
    const stk::mesh::Entity slaveNode = vecPair.second;
    // pointer to data
    double *masterField = (double *)stk::mesh::field_data(*theField, masterNode);
    double *slaveField = (double *)stk::mesh::field_data(*theField, slaveNode);

    for ( unsigned j = 0; j < sizeOfField; ++j ) {
      const double maxValue = std::max(masterField[j],slaveField[j]);
      masterField[j] = maxValue; 
      slaveField[j] = maxValue;
    }
  }

  // parallel communicate shared and aura-ed entities
  parallel_communicate_field(theField);

}

//--------------------------------------------------------------------------
//-------- add_slave_to_master ---------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::add_slave_to_master(
  stk::mesh::FieldBase *theField,
  const unsigned &sizeOfField,
  const bool &bypassFieldCheck)
{
  // iterate vector of masterEntity:slaveEntity pairs
  if ( bypassFieldCheck ) {
    // fields are expected to be defined on all master/slave nodes
    for ( size_t k = 0; k < masterSlaveCommunicator_.size(); ++k) {
      // extract master node and slave node
      EntityPair vecPair = masterSlaveCommunicator_[k];
      const stk::mesh::Entity masterNode = vecPair.first;
      const stk::mesh::Entity slaveNode = vecPair.second;
      // pointer to data
      double *masterField = (double *)stk::mesh::field_data(*theField, masterNode);
      const double *slaveField = (double *)stk::mesh::field_data(*theField, slaveNode);
      // add in contribution
      for ( unsigned j = 0; j < sizeOfField; ++j ) {
	masterField[j] += slaveField[j];
      }
    }
  }
  else {
    // more costly check to see if fields are defined on master/slave nodes    
    for ( size_t k = 0; k < masterSlaveCommunicator_.size(); ++k) {      
      // extract master node and slave node
      EntityPair vecPair = masterSlaveCommunicator_[k];
      const stk::mesh::Entity masterNode = vecPair.first;
      const stk::mesh::Entity slaveNode = vecPair.second;
      // pointer to data
      double *masterField = (double *)stk::mesh::field_data(*theField, masterNode);
      if ( NULL != masterField ) {
	const double *slaveField = (double *)stk::mesh::field_data(*theField, slaveNode);	
	// add in contribution
	for ( unsigned j = 0; j < sizeOfField; ++j ) {
	  masterField[j] += slaveField[j];
	}
      }
    }
  }

}

//--------------------------------------------------------------------------
//-------- set_slave_to_master ---------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::set_slave_to_master(
  stk::mesh::FieldBase *theField,
  const unsigned &sizeOfField,
  const bool &bypassFieldCheck)
{
  // iterate vector of masterEntity:slaveEntity pairs
  if ( bypassFieldCheck ) {
    // fields are expected to be defined on all master/slave nodes
    for ( size_t k = 0; k < masterSlaveCommunicator_.size(); ++k) {
      // extract master node and slave node
      EntityPair vecPair = masterSlaveCommunicator_[k];
      const stk::mesh::Entity masterNode = vecPair.first;
      const stk::mesh::Entity slaveNode = vecPair.second;
      // pointer to data
      const double *masterField = (double *)stk::mesh::field_data(*theField, masterNode);
      double *slaveField = (double *)stk::mesh::field_data(*theField, slaveNode);
      // set master to slave
      for ( unsigned j = 0; j < sizeOfField; ++j ) {
	slaveField[j] = masterField[j];
      }
    }
  }
  else {
    // more costly check to see if fields are defined on master/slave nodes    
    for ( size_t k = 0; k < masterSlaveCommunicator_.size(); ++k) {
      // extract master node and slave node
      EntityPair vecPair = masterSlaveCommunicator_[k];
      const stk::mesh::Entity masterNode = vecPair.first;
      const stk::mesh::Entity slaveNode = vecPair.second;
      // pointer to data
      const double *masterField = (double *)stk::mesh::field_data(*theField, masterNode);
      
      if ( NULL != masterField ) {
	double *slaveField = (double *)stk::mesh::field_data(*theField, slaveNode);
	// set master to slave
	for ( unsigned j = 0; j < sizeOfField; ++j ) {
	  slaveField[j] = masterField[j];
	}
      }
    }
  }
}

} // namespace nalu
} // namespace sierra
