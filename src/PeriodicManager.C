/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <PeriodicManager.h>
#include <NaluEnv.h>
#include <Realm.h>
#include <utils/StkHelpers.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SortAndUnique.hpp>

// stk_search
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>

// vector
#include <vector>
#include <map>
#include <string>

namespace sierra{
namespace nalu{

PeriodicManager::PeriodicManager(
   Realm &realm)
  : realm_(realm ),
    searchTolerance_(1.0e-8),
    periodicGhosting_(NULL),
    ghostingName_("nalu_periodic"),
    timerSearch_(0.0),
    errorCount_(0),
    maxErrorCount_(10),
    amplificationFactor_(5.0)
{
  // does nothing
}

PeriodicManager::~PeriodicManager()
{
  // nothing to delete
}

//--------------------------------------------------------------------------
//-------- initialize_error_count ------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::initialize_error_count()
{
  errorCount_ = 0;
}

//--------------------------------------------------------------------------
//-------- add periodic pair -----------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::add_periodic_pair(
  stk::mesh::Part * monarchMeshPart,
  stk::mesh::Part * subjectMeshPart,
  const double &userSearchTolerance,
  const std::string &searchMethodName)
{
  // use most stringent tolerance (min) for all of user specifications
  searchTolerance_ = std::min(searchTolerance_, userSearchTolerance);

  // form the subject part vector
  subjectPartVector_.push_back(subjectMeshPart);

  // form the selector pair
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::Selector monarchSelect = meta_data.locally_owned_part() & stk::mesh::Selector(*monarchMeshPart);
  stk::mesh::Selector subjectSelect = meta_data.locally_owned_part() & stk::mesh::Selector(*subjectMeshPart);
  SelectorPair periodicSelectorPair(monarchSelect, subjectSelect);
  periodicSelectorPairs_.push_back(periodicSelectorPair);

  // determine search method for this pair; default is stk_kdtree
  stk::search::SearchMethod searchMethod = stk::search::KDTREE;
  if ( searchMethodName != "stk_kdtree" )
    NaluEnv::self().naluOutputP0() << "PeriodicManager::search_method only supports stk_kdtree" 
                                   << std::endl;
  searchMethodVec_.push_back(searchMethod);
}

//--------------------------------------------------------------------------
//-------- get_subject_part_vector ----------------------------------------------
//--------------------------------------------------------------------------
const stk::mesh::PartVector &
PeriodicManager::get_subject_part_vector()
{
  return subjectPartVector_;
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
  if ( periodicSelectorPairs_.size() == 0 )
     throw std::runtime_error("PeriodiocBC::Error: No periodic pair provided");

  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "Periodic Review:  realm: " << realm_.name_ << std::endl;
  NaluEnv::self().naluOutputP0() << "=========================" << std::endl;

  augment_periodic_selector_pairs();

  // initialize translation and rotation vectors
  initialize_translation_vector();

  // translate
  for ( size_t k = 0; k < periodicSelectorPairs_.size(); ++k) {
    determine_translation(periodicSelectorPairs_[k].first, periodicSelectorPairs_[k].second,
        translationVector_[k], rotationVector_[k]);
  }

  remove_redundant_subject_nodes();

  // search and constraint mapping
  finalize_search();

  // provide Nalu id update
  update_global_id_field();
}

//--------------------------------------------------------------------------
//-------- augment_periodic_selector_pairs ---------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::augment_periodic_selector_pairs()
{
  
  const size_t pairSize = periodicSelectorPairs_.size();

  switch ( pairSize ) {

    case 1:
      break; // nothing to do

    case 2: {

      // monarch/subject selectors
      stk::mesh::Selector &monarchA = periodicSelectorPairs_[0].first;
      stk::mesh::Selector &monarchB = periodicSelectorPairs_[1].first;

      stk::mesh::Selector &subjectA = periodicSelectorPairs_[0].second;
      stk::mesh::Selector &subjectB = periodicSelectorPairs_[1].second;

      // push back intersection selector pairs
      periodicSelectorPairs_.push_back(std::make_pair(monarchA & monarchB, subjectA & subjectB));

      // need a search method; arbitrarily choose the first method specified
      stk::search::SearchMethod searchMethod = searchMethodVec_[0];
      searchMethodVec_.push_back(searchMethod);
    
      break;
    }

    case 3: {
      const stk::mesh::Selector monarchA = periodicSelectorPairs_[0].first;
      const stk::mesh::Selector monarchB = periodicSelectorPairs_[1].first;
      const stk::mesh::Selector monarchC = periodicSelectorPairs_[2].first;

      const stk::mesh::Selector subjectA = periodicSelectorPairs_[0].second;
      const stk::mesh::Selector subjectB = periodicSelectorPairs_[1].second;
      const stk::mesh::Selector subjectC = periodicSelectorPairs_[2].second;

      // push back intersection selector pairs
      periodicSelectorPairs_.push_back(std::make_pair(monarchA & monarchB, subjectA & subjectB));
      periodicSelectorPairs_.push_back(std::make_pair(monarchB & monarchC, subjectB & subjectC));
      periodicSelectorPairs_.push_back(std::make_pair(monarchA & monarchC, subjectA & subjectC));
      periodicSelectorPairs_.push_back(std::make_pair(monarchA & monarchB & monarchC, subjectA & subjectB & subjectC));

      // need a search method; arbitrarily choose the first method specified
      stk::search::SearchMethod searchMethod = searchMethodVec_[0];
      searchMethodVec_.push_back(searchMethod); // 3
      searchMethodVec_.push_back(searchMethod); // 4
      searchMethodVec_.push_back(searchMethod); // 5
      searchMethodVec_.push_back(searchMethod); // 6

      break;
    }

    default: {
      NaluEnv::self().naluOutputP0() << "more than three periodic pairs assumes no common subject nodes " << std::endl;
      break;
    }
  }
}

//--------------------------------------------------------------------------
//-------- initialize_translation_vector -----------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::initialize_translation_vector()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  translationVector_.resize(periodicSelectorPairs_.size());
  rotationVector_.resize(periodicSelectorPairs_.size());
  const int nDim = meta_data.spatial_dimension();
  for ( size_t k = 0; k < periodicSelectorPairs_.size(); ++k ) {
    translationVector_[k].resize(nDim);
    rotationVector_[k].resize(nDim);
  }
}

//--------------------------------------------------------------------------
//-------- determine_translation -------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::determine_translation(
    stk::mesh::Selector monarchSelector,
    stk::mesh::Selector subjectSelector,
    std::vector<double> &translationVector,
    std::vector<double> &rotationVector)
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // fields
  VectorFieldType *coordinates = meta_data.get_field<double>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  const int nDim = meta_data.spatial_dimension();

  // Monarch: global_sum_coords_monarch
  std::vector<double> local_sum_coords_monarch(nDim, 0.0), global_sum_coords_monarch(nDim, 0.0);
  size_t numberMonarchNodes = 0, g_numberMonarchNodes = 0;

  stk::mesh::BucketVector const& monarch_node_buckets = realm_.get_buckets( stk::topology::NODE_RANK, monarchSelector);

  for ( stk::mesh::BucketVector::const_iterator ib = monarch_node_buckets.begin();
        ib != monarch_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length   = b.size();
    // point to data
    const double * coords = stk::mesh::field_data(*coordinates, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      numberMonarchNodes += 1;

      // define offset for all nodal fields that are of nDim
      const size_t offSet = k*nDim;

      // sum local coords for translation; define localCoords for bounding point
      for (int j = 0; j < nDim; ++j ) {
        const double cxj = coords[offSet+j];
        local_sum_coords_monarch[j] += cxj;
      }
    }
  }
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &local_sum_coords_monarch[0], &global_sum_coords_monarch[0], nDim);
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &numberMonarchNodes, &g_numberMonarchNodes, 1);

  // Subject: global_sum_coords_subject
  std::vector<double> local_sum_coords_subject(nDim, 0.0), global_sum_coords_subject(nDim, 0.0);
  size_t numberSubjectNodes = 0, g_numberSubjectNodes = 0;

  stk::mesh::BucketVector const& subject_node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, subjectSelector);

  for ( stk::mesh::BucketVector::const_iterator ib = subject_node_buckets.begin();
        ib != subject_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length   = b.size();
    // point to data
    const double * coords = stk::mesh::field_data(*coordinates, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      numberSubjectNodes += 1;
      // define offset for all nodal fields that are of nDim
      const size_t offSet = k*nDim;
      for (int j = 0; j < nDim; ++j ) {
        local_sum_coords_subject[j] += coords[offSet+j];
      }
    }
  }
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &local_sum_coords_subject[0], &global_sum_coords_subject[0], nDim);
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &numberSubjectNodes, &g_numberSubjectNodes, 1);

  // save off translation and rotation
  for (int j = 0; j < nDim; ++j ) {
    translationVector[j] = (global_sum_coords_monarch[j] - global_sum_coords_subject[j]) / g_numberMonarchNodes;
    rotationVector[j] = global_sum_coords_monarch[j] / g_numberMonarchNodes;
  }

  NaluEnv::self().naluOutputP0() << "Translating [ ";
  for (int j = 0; j < nDim; ++j ) {  NaluEnv::self().naluOutputP0() << translationVector[j] << " "; }
  NaluEnv::self().naluOutputP0() << "] Monarch/Subject pair " << std::endl;
}

//--------------------------------------------------------------------------
//-------- remove_redundant_subject_nodes ----------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::remove_redundant_subject_nodes()
{

  switch (periodicSelectorPairs_.size()) {

    case 1: case 2:
      break; // nothing to do

    case 3: {
      // subject selectors
      stk::mesh::Selector &subjectA = periodicSelectorPairs_[0].second;
      stk::mesh::Selector &subjectB = periodicSelectorPairs_[1].second;

      // intersection of A/B
      stk::mesh::Selector subjectIntersection = subjectA & subjectB;

      // now remove redundant [corner] nodes
      periodicSelectorPairs_[0].second = subjectA - subjectIntersection;
      periodicSelectorPairs_[1].second = subjectB - subjectIntersection;

      break;
    }

    case 7: {
      // subject selectors
      const stk::mesh::Selector subjectA = periodicSelectorPairs_[0].second;
      const stk::mesh::Selector subjectB = periodicSelectorPairs_[1].second;
      const stk::mesh::Selector subjectC = periodicSelectorPairs_[2].second;

      // intersection of A/B/C (corner nodes)
      const stk::mesh::Selector subjectABC = subjectA & subjectB & subjectC;

      // intersection of A/B/C (edges of box)
      const stk::mesh::Selector subjectAB = subjectA & subjectB;
      const stk::mesh::Selector subjectAC = subjectA & subjectC;
      const stk::mesh::Selector subjectBC = subjectB & subjectC;

      // now remove redundant [corner] nodes
      periodicSelectorPairs_[0].second = subjectA - (subjectAB | subjectAC);
      periodicSelectorPairs_[1].second = subjectB - (subjectAB | subjectBC);
      periodicSelectorPairs_[2].second = subjectC - (subjectAC | subjectBC);

      // now remove redundant [edges of box] nodes
      periodicSelectorPairs_[3].second = subjectAB - subjectABC;
      periodicSelectorPairs_[4].second = subjectBC - subjectABC;
      periodicSelectorPairs_[5].second = subjectAC - subjectABC;

      break;
    }
    default: {
      NaluEnv::self().naluOutputP0() << "more than three periodic pairs assumes no common subject nodes " << std::endl;
      break;
    }
  }
}

//--------------------------------------------------------------------------
//-------- finalize_search -------------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::finalize_search()
{
  // clear vectors
  searchKeyVector_.clear();
  monarchSubjectCommunicator_.clear();

  // process each pair
  for ( size_t k = 0; k < periodicSelectorPairs_.size(); ++k) {
    populate_search_key_vec(periodicSelectorPairs_[k].first, periodicSelectorPairs_[k].second,
                            translationVector_[k], searchMethodVec_[k]);
  }
  
  // manage ghosting
  manage_ghosting_object();

  // check to see if we need to attempt once more
  error_check();
}

//--------------------------------------------------------------------------
//-------- populate_search_key_vec -----------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::populate_search_key_vec(
    stk::mesh::Selector monarchSelector,
    stk::mesh::Selector subjectSelector,
    std::vector<double> &translationVector,
    const stk::search::SearchMethod searchMethod)
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  // required data structures; monarch/subject
  std::vector<sphereBoundingBox> sphereBoundingBoxMonarchVec;
  std::vector<sphereBoundingBox> sphereBoundingBoxSubjectVec;

  // fields
  VectorFieldType *coordinates = meta_data.get_field<double>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  const int nDim = meta_data.spatial_dimension();

  // Point
  const double pointRadius = searchTolerance_;
  Point monarchCenter, subjectCenter;

  // Monarch: setup sphereBoundingBoxMonarchVec,
  stk::mesh::BucketVector const& monarch_node_buckets = realm_.get_buckets( stk::topology::NODE_RANK, monarchSelector);

  for ( stk::mesh::BucketVector::const_iterator ib = monarch_node_buckets.begin();
        ib != monarch_node_buckets.end() ; ++ib ) {
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

      // sum local coords for translation; define localCoords for bounding point
      for (int j = 0; j < nDim; ++j ) {
        const double cxj = coords[offSet+j];
        monarchCenter[j] = cxj;
      }

      // create the bounding point sphere and push back
      sphereBoundingBox theSphere( Sphere(monarchCenter, pointRadius), theIdent);
      sphereBoundingBoxMonarchVec.push_back(theSphere);
    }
  }

  // Subject: setup sphereBoundingBoxSubjectVec; translate subject onto monarch
  stk::mesh::BucketVector const& subject_node_buckets =
  realm_.get_buckets( stk::topology::NODE_RANK, subjectSelector);
  for ( stk::mesh::BucketVector::const_iterator ib = subject_node_buckets.begin();
       ib != subject_node_buckets.end() ; ++ib ) {
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
        subjectCenter[j] = xj + translationVector[j];
      }
      sphereBoundingBox theSphere(Sphere(subjectCenter, pointRadius), theIdent);
      sphereBoundingBoxSubjectVec.push_back(theSphere);
    }
  }

  // will want to stuff product of search to a single vector
  std::vector<std::pair<theEntityKey, theEntityKey> > searchKeyPair;
  double timeA = NaluEnv::self().nalu_time();
  stk::search::coarse_search(sphereBoundingBoxSubjectVec, sphereBoundingBoxMonarchVec, searchMethod, NaluEnv::self().parallel_comm(), searchKeyPair);
  timerSearch_ += (NaluEnv::self().nalu_time() - timeA);

  // populate searchKeyVector_; culmination of all monarch/subjects
  searchKeyVector_.insert(searchKeyVector_.end(), searchKeyPair.begin(), searchKeyPair.end());
}

//--------------------------------------------------------------------------
//-------- error_check -----------------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::error_check()
{
  // number of subject nodes should equal the size of the searchKeyVector_
  size_t l_totalNumber[2] = {0,0};

  // extract total locally owned subject nodes
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::Selector s_locally_owned = meta_data.locally_owned_part()
    & stk::mesh::selectUnion(subjectPartVector_);
  stk::mesh::BucketVector const& subject_node_buckets = realm_.get_buckets( stk::topology::NODE_RANK, s_locally_owned);
  for ( stk::mesh::BucketVector::const_iterator ib = subject_node_buckets.begin();
	ib != subject_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length   = b.size();
    l_totalNumber[0] += length;
  }
  
  // extract locally owned subject nodes from the search
  for (size_t i=0, size=searchKeyVector_.size(); i<size; ++i) {
    if ( NaluEnv::self().parallel_rank() == searchKeyVector_[i].second.proc())
      l_totalNumber[1] += 1;
  }
  
  // parallel sum and check
  size_t g_totalNumber[2] = {0,0};
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), l_totalNumber, g_totalNumber, 2);

  // hard error check
  if ( g_totalNumber[0] != g_totalNumber[1]) {
    // increment and report
    errorCount_++;
    NaluEnv::self().naluOutputP0() << "_____________________________________" << std::endl;
    NaluEnv::self().naluOutputP0() << "Probable issue with Search on attempt: " << errorCount_ << std::endl;
    NaluEnv::self().naluOutputP0() << "the total number of subject nodes (" << g_totalNumber[0] << ")" << std::endl;
    NaluEnv::self().naluOutputP0() << "does not equal the product of the search (" << g_totalNumber[1] << ")" << std::endl;

    // check to see if we should try again..
    if ( errorCount_ == maxErrorCount_ ) {
      NaluEnv::self().naluOutputP0() << "ABORT: Too many attempts; please check your mesh" << std::endl;
      throw std::runtime_error("PeriodiocBC::Error: Too many attempts; please check your mesh");
    }
    // reduce or increase search tolerance based on number of subjects and total
    if ( g_totalNumber[0] > g_totalNumber[1] ) {
      searchTolerance_ *= amplificationFactor_; // subject > total; increase tolerance
      NaluEnv::self().naluOutputP0() << "The algorithm will increase the search tolerance: " << searchTolerance_ << std::endl;
    }
    else {
      searchTolerance_ /= amplificationFactor_; // subject < total; reduce tolerance
      NaluEnv::self().naluOutputP0() << "The algorithm will reduce the search tolerance: " << searchTolerance_ << std::endl;
    }

    // try again
    finalize_search();
  }
  else {
    NaluEnv::self().naluOutputP0() << "---------------------------------------------------" << std::endl;
    NaluEnv::self().naluOutputP0() << "Parallel consistency noted in monarch/subject pairings: "
                                   << g_totalNumber[0] << "/"<< g_totalNumber[1] << std::endl;
    NaluEnv::self().naluOutputP0() << "---------------------------------------------------" << std::endl;
    NaluEnv::self().naluOutputP0() << std::endl;
  }
}

//--------------------------------------------------------------------------
//-------- manage_ghosting_object ------------------------------------------
//--------------------------------------------------------------------------
void add_range_nodes_to_sharers_of_domain_nodes(stk::mesh::BulkData& bulk_data,
                                                  const PeriodicManager::SearchKeyVector& searchKeyVector,
                                                  std::vector<stk::mesh::EntityProc>& sendNodes)
{
    stk::CommSparse commSparse(bulk_data.parallel());

    auto packingLambda = [&]() {
        int theRank = NaluEnv::self().parallel_rank();
        std::vector<int> sharingProcs;
        for(size_t i=0; i<searchKeyVector.size(); ++i) {
            int domainProc = searchKeyVector[i].first.proc();
            int rangeProc = searchKeyVector[i].second.proc();

            stk::mesh::EntityKey domainKey = searchKeyVector[i].first.id();
            stk::mesh::Entity domainNode = bulk_data.get_entity(domainKey);

            if ((theRank == domainProc) && bulk_data.bucket(domainNode).shared()) {
                stk::mesh::EntityId rangeId = searchKeyVector[i].second.id().id();
                stk::CommBuffer& sbuf = commSparse.send_buffer(rangeProc);

                bulk_data.comm_shared_procs(domainKey, sharingProcs);
                if (theRank == rangeProc) {
                    stk::mesh::Entity rangeNode = bulk_data.get_entity(stk::topology::NODE_RANK, rangeId);
                    for(int p : sharingProcs) {
                        sendNodes.push_back(stk::mesh::EntityProc(rangeNode, p));
                    }
                }
                else {
                    for(int p : sharingProcs) {
                        if (p != theRank && p != rangeProc) {
                            sbuf.pack(rangeId);
                            sbuf.pack(p);
                        }
                    }
                }
            }
        }
    };

    stk::pack_and_communicate(commSparse, packingLambda);

    stk::unpack_communications(commSparse, [&](int p)
    {
        stk::CommBuffer& rbuf = commSparse.recv_buffer(p);
        stk::mesh::EntityId gid;
        rbuf.unpack(gid);
        int proc;
        rbuf.unpack(proc);
        sendNodes.push_back(stk::mesh::EntityProc(bulk_data.get_entity(stk::topology::NODE_RANK, gid), proc));
    });
}

void
PeriodicManager::manage_ghosting_object()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  unsigned theRank = NaluEnv::self().parallel_rank();

  std::vector<stk::mesh::EntityProc> sendNodes;
  std::vector<stk::mesh::EntityProc> sendElems;
  std::vector<int> sharingProcs;
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

  add_range_nodes_to_sharers_of_domain_nodes(bulk_data, searchKeyVector_, sendNodes);

  size_t numNodes = sendNodes.size();
  size_t g_numNodes = 0;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &numNodes, &g_numNodes, 1);
  if ( g_numNodes > 0) {
    // check if we need to ghost
    bulk_data.modification_begin();
    if ( periodicGhosting_ == NULL )
      periodicGhosting_ = &bulk_data.create_ghosting(ghostingName_);
    else
      bulk_data.destroy_ghosting(*periodicGhosting_);
    bulk_data.change_ghosting(*periodicGhosting_, sendNodes);
    bulk_data.modification_end();

    populate_ghost_comm_procs(bulk_data, *periodicGhosting_, ghostCommProcs_);
  }

  // now populate monarch subject communicator
  for (size_t i=0, size=searchKeyVector_.size(); i<size; ++i) {
     stk::mesh::Entity domainNode = bulk_data.get_entity(searchKeyVector_[i].first.id());
     stk::mesh::Entity rangeNode = bulk_data.get_entity(searchKeyVector_[i].second.id());
     // unique monarch:subject communicator
     EntityPair theFirstPair = std::make_pair(rangeNode, domainNode);
     monarchSubjectCommunicator_.push_back(theFirstPair);
  }
}

//--------------------------------------------------------------------------
//-------- get_ghosting_object ---------------------------------------------
//--------------------------------------------------------------------------
stk::mesh::Ghosting *
PeriodicManager::get_ghosting_object()
{
  return periodicGhosting_;
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
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
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

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  // no need to update periodically ghosted fields..

  // vector of monarchEntity:subjectEntity pairs
  for ( size_t k = 0; k < monarchSubjectCommunicator_.size(); ++k) {

    // extract monarch node and subject node
    EntityPair vecPair = monarchSubjectCommunicator_[k];
    const stk::mesh::Entity monarchNode = vecPair.first;
    const stk::mesh::Entity subjectNode = vecPair.second;

    // pointer to data
    const stk::mesh::EntityId monarchGlobalId = bulk_data.identifier(monarchNode );
    stk::mesh::EntityId *subjectGlobalId = stk::mesh::field_data(*realm_.naluGlobalId_, subjectNode);

    // set new value
    *subjectGlobalId = monarchGlobalId;

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
  const unsigned sizeOfField,
  const bool bypassFieldCheck,
  const bool addSubjects,
  const bool setSubjects)
{
  // update periodically ghosted fields within add_ and set_
  if ( addSubjects )
    add_subject_to_monarch(theField, sizeOfField, bypassFieldCheck);
  if ( setSubjects )
    set_subject_to_monarch(theField, sizeOfField, bypassFieldCheck);

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

  periodic_parallel_communicate_field(theField);

  for ( size_t k = 0; k < monarchSubjectCommunicator_.size(); ++k) {
    // extract monarch node and subject node
    EntityPair vecPair = monarchSubjectCommunicator_[k];
    const stk::mesh::Entity monarchNode = vecPair.first;
    const stk::mesh::Entity subjectNode = vecPair.second;
    // pointer to data
    double *monarchField = (double *)stk::mesh::field_data(*theField, monarchNode);
    double *subjectField = (double *)stk::mesh::field_data(*theField, subjectNode);

    for ( unsigned j = 0; j < sizeOfField; ++j ) {
      const double maxValue = std::max(monarchField[j],subjectField[j]);
      monarchField[j] = maxValue; 
      subjectField[j] = maxValue;
    }
  }

  // parallel communicate shared and aura-ed entities
  parallel_communicate_field(theField);

}

//--------------------------------------------------------------------------
//-------- add_subject_to_monarch ------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::add_subject_to_monarch(
  stk::mesh::FieldBase *theField,
  const unsigned &sizeOfField,
  const bool &bypassFieldCheck)
{
  
  periodic_parallel_communicate_field(theField);

  // iterate vector of monarchEntity:subjectEntity pairs
  if ( bypassFieldCheck ) {
    // fields are expected to be defined on all monarch/subject nodes
    for ( size_t k = 0; k < monarchSubjectCommunicator_.size(); ++k) {
      // extract monarch node and subject node
      EntityPair vecPair = monarchSubjectCommunicator_[k];
      const stk::mesh::Entity monarchNode = vecPair.first;
      const stk::mesh::Entity subjectNode = vecPair.second;
      // pointer to data
      double *monarchField = (double *)stk::mesh::field_data(*theField, monarchNode);
      const double *subjectField = (double *)stk::mesh::field_data(*theField, subjectNode);
      // add in contribution
      for ( unsigned j = 0; j < sizeOfField; ++j ) {
        monarchField[j] += subjectField[j];
      }
    }
  }
  else {
    // more costly check to see if fields are defined on monarch/subject nodes    
    for ( size_t k = 0; k < monarchSubjectCommunicator_.size(); ++k) {      
      // extract monarch node and subject node
      EntityPair vecPair = monarchSubjectCommunicator_[k];
      const stk::mesh::Entity monarchNode = vecPair.first;
      const stk::mesh::Entity subjectNode = vecPair.second;
      // pointer to data
      double *monarchField = (double *)stk::mesh::field_data(*theField, monarchNode);
      if ( NULL != monarchField ) {
        const double *subjectField = (double *)stk::mesh::field_data(*theField, subjectNode);
        // add in contribution
        for ( unsigned j = 0; j < sizeOfField; ++j ) {
          monarchField[j] += subjectField[j];
        }
      }
    }
  }

  periodic_parallel_communicate_field(theField);

}

//--------------------------------------------------------------------------
//-------- set_subject_to_monarch ------------------------------------------
//--------------------------------------------------------------------------
void
PeriodicManager::set_subject_to_monarch(
  stk::mesh::FieldBase *theField,
  const unsigned &sizeOfField,
  const bool &bypassFieldCheck)
{

  periodic_parallel_communicate_field(theField);

  // iterate vector of monarchEntity:subjectEntity pairs
  if ( bypassFieldCheck ) {
    // fields are expected to be defined on all monarch/subject nodes
    for ( size_t k = 0; k < monarchSubjectCommunicator_.size(); ++k) {
      // extract monarch node and subject node
      EntityPair vecPair = monarchSubjectCommunicator_[k];
      const stk::mesh::Entity monarchNode = vecPair.first;
      const stk::mesh::Entity subjectNode = vecPair.second;
      // pointer to data
      const double *monarchField = (double *)stk::mesh::field_data(*theField, monarchNode);
      double *subjectField = (double *)stk::mesh::field_data(*theField, subjectNode);
      // set monarch to subject
      for ( unsigned j = 0; j < sizeOfField; ++j ) {
        subjectField[j] = monarchField[j];
      }
    }
  }
  else {
    // more costly check to see if fields are defined on monarch/subject nodes    
    for ( size_t k = 0; k < monarchSubjectCommunicator_.size(); ++k) {
      // extract monarch node and subject node
      EntityPair vecPair = monarchSubjectCommunicator_[k];
      const stk::mesh::Entity monarchNode = vecPair.first;
      const stk::mesh::Entity subjectNode = vecPair.second;
      // pointer to data
      const double *monarchField = (double *)stk::mesh::field_data(*theField, monarchNode);
      
      if ( NULL != monarchField ) {
        double *subjectField = (double *)stk::mesh::field_data(*theField, subjectNode);
        // set monarch to subject
        for ( unsigned j = 0; j < sizeOfField; ++j ) {
          subjectField[j] = monarchField[j];
        }
      }
    }
  }

  periodic_parallel_communicate_field(theField);

}

} // namespace nalu
} // namespace sierra
