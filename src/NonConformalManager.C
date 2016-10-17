/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <NonConformalInfo.h>
#include <NonConformalManager.h>
#include <master_element/MasterElement.h>
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

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/util/SortAndUnique.hpp>

// vector and pair
#include <vector>

namespace sierra{
namespace nalu{

class HaloInfo;

//==========================================================================
// Class Definition
//==========================================================================
// Acon_NonConformalManager - manages nonConformal info
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
NonConformalManager::NonConformalManager(
  Realm &realm,
  const bool ncAlgDetailedOutput)
  : realm_(realm ),
    ncAlgDetailedOutput_(ncAlgDetailedOutput),
    nonConformalGhosting_(NULL)
{
  // do nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
NonConformalManager::~NonConformalManager()
{
  // delete nonConformal info objects:
  std::vector<NonConformalInfo *>::iterator ic;
  for( ic=nonConformalInfoVec_.begin();
       ic!=nonConformalInfoVec_.end(); ++ic )
    delete (*ic);
}

void combineElemsAndDownwardRelations(const stk::mesh::BulkData& bulk,
                                  const stk::mesh::EntityProcVec& elemsToGhost,
                                  stk::mesh::EntityProcVec& newSendGhosts)
{
  newSendGhosts.reserve(elemsToGhost.size());
  for(const stk::mesh::EntityProc& elemProc : elemsToGhost) {
    newSendGhosts.push_back(elemProc);
    stk::mesh::EntityRank thisRank = bulk.entity_rank(elemProc.first);
    for(stk::mesh::EntityRank irank=stk::topology::NODE_RANK; irank<thisRank; ++irank) {
      unsigned num = bulk.num_connectivity(elemProc.first, irank);
      const stk::mesh::Entity* entities = bulk.begin(elemProc.first, irank);
      for(unsigned i=0; i<num; ++i) {
        newSendGhosts.push_back(stk::mesh::EntityProc(entities[i], elemProc.second));
      }
    }
  }
  stk::util::sort_and_unique(newSendGhosts);
}

void keepElemsNotAlreadyGhosted(const stk::mesh::BulkData& bulk,
                                const stk::mesh::EntityProcVec& newSendGhosts,
                                const stk::mesh::EntityProcVec& intersection,
                                stk::mesh::EntityProcVec& elemsToGhost)
{
  elemsToGhost.clear();
  for(const stk::mesh::EntityProc& entityProc : newSendGhosts) {
    if (bulk.entity_rank(entityProc.first) == stk::topology::ELEM_RANK) {
      if (!std::binary_search(intersection.begin(), intersection.end(), entityProc)) {
        elemsToGhost.push_back(entityProc);
      }
    }
  }
}

void fillSendGhostsToRemoveFromGhosting(const stk::mesh::EntityProcVec& curSendGhosts,
                                        const stk::mesh::EntityProcVec& intersection,
                                        stk::mesh::EntityProcVec& sendGhostsToRemove)
{
  sendGhostsToRemove.reserve(curSendGhosts.size() - intersection.size());
  for(size_t i=0; i<curSendGhosts.size(); ++i) {
    if (!std::binary_search(intersection.begin(), intersection.end(), curSendGhosts[i])) {
      sendGhostsToRemove.push_back(curSendGhosts[i]);
    }
  }
}

void communicateToFillRecvGhostsToRemove(const stk::mesh::BulkData& bulk,
                                         const stk::mesh::EntityProcVec& sendGhostsToRemove,
                                         std::vector<stk::mesh::EntityKey>& recvGhostsToRemove)
{
  stk::CommSparse commSparse(bulk.parallel());
  stk::pack_and_communicate(commSparse, [&]() {
    for(const stk::mesh::EntityProc& entityProc : sendGhostsToRemove) {
      stk::mesh::EntityKey key = bulk.entity_key(entityProc.first);
      stk::CommBuffer& buf = commSparse.send_buffer(entityProc.second);
      buf.pack<stk::mesh::EntityKey>(key);
    }
  });

  int numProcs = bulk.parallel_size();
  for(int p=0; p<numProcs; ++p) {
    if (p == bulk.parallel_rank()) {
      continue;
    }
    stk::CommBuffer& buf = commSparse.recv_buffer(p);
    while(buf.remaining()) {
      stk::mesh::EntityKey key;
      buf.unpack<stk::mesh::EntityKey>(key);
      recvGhostsToRemove.push_back(key);
    }
  }
}

void
NonConformalManager::computePreciseGhostingLists(const stk::mesh::BulkData& bulk,
                                                 stk::mesh::EntityProcVec& elemsToGhost,
                                                 stk::mesh::EntityProcVec& curSendGhosts,
                                                 std::vector<stk::mesh::EntityKey>& recvGhostsToRemove)
{
  stk::util::sort_and_unique(curSendGhosts);
  stk::util::sort_and_unique(elemsToGhost);

  stk::mesh::EntityProcVec newSendGhosts;
  combineElemsAndDownwardRelations(bulk, elemsToGhost, newSendGhosts);

  stk::mesh::EntityProcVec intersection;
  std::set_intersection(curSendGhosts.begin(), curSendGhosts.end(),
                        newSendGhosts.begin(), newSendGhosts.end(),
                        std::back_inserter(intersection));

  keepElemsNotAlreadyGhosted(bulk, newSendGhosts, intersection, elemsToGhost);

  stk::mesh::EntityProcVec sendGhostsToRemove;
  fillSendGhostsToRemoveFromGhosting(curSendGhosts, intersection, sendGhostsToRemove);

  communicateToFillRecvGhostsToRemove(bulk, sendGhostsToRemove, recvGhostsToRemove);
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
NonConformalManager::initialize()
{

  const double timeA = NaluEnv::self().nalu_time();

  elemsToGhost_.clear();

  std::vector<stk::mesh::EntityKey> recvGhostsToRemove;
  stk::mesh::EntityProcVec currentSendGhosts;

  if (nonConformalGhosting_ != NULL) {
      nonConformalGhosting_->send_list(currentSendGhosts);
  }
  
  // loop over nonConformalInfo and initialize to update the elemsToGhost_ vector.
  for ( size_t k = 0; k < nonConformalInfoVec_.size(); ++k )
    nonConformalInfoVec_[k]->initialize();
 
  //We want elemsToGhost_ to only contain elements not already ghosted, and
  //we want a list of receive-ghosts that no longer need to be ghosted.
  computePreciseGhostingLists(realm_.bulk_data(), elemsToGhost_, currentSendGhosts, recvGhostsToRemove);

  // check for ghosting need
  size_t local[2] = {elemsToGhost_.size(), recvGhostsToRemove.size()};
  size_t global[2] = {0, 0};
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), local, global, 2);

  if (global[0] > 0 || global[1] > 0) {
      NaluEnv::self().naluOutputP0() << "NonConformal alg will ghost a number of entities: "
                    << global[0]<<" and remove "<<global[1]<<" entities from ghosting."  << std::endl;

      manage_ghosting(recvGhostsToRemove);
  }
  else {
    NaluEnv::self().naluOutputP0() << "NonConformal alg will NOT ghost entities: " << std::endl;
  }

  // complete search
  for ( size_t k = 0; k < nonConformalInfoVec_.size(); ++k )
    nonConformalInfoVec_[k]->complete_search();

  // provide diagnosis
  if ( ncAlgDetailedOutput_ ) {
    for ( size_t k = 0; k < nonConformalInfoVec_.size(); ++k )
      nonConformalInfoVec_[k]->provide_diagnosis();
  }

  if (nonConformalGhosting_ != nullptr) {
    //refresh all fields on ghosts. potential optimization: only refresh needed fields?
    const stk::mesh::FieldVector& fields = realm_.bulk_data().mesh_meta_data().get_fields();
    std::vector<const stk::mesh::FieldBase*> const_fields;
    const_fields.reserve(fields.size());
    for(stk::mesh::FieldBase* fieldptr : fields) {
      const_fields.push_back(fieldptr);
    }
    stk::mesh::communicate_field_data(*nonConformalGhosting_, const_fields);
  }

  // end time
  const double timeB = NaluEnv::self().nalu_time();
  realm_.timerContact_ += (timeB-timeA);
}

//--------------------------------------------------------------------------
//-------- manage_ghosting -------------------------------------------------
//--------------------------------------------------------------------------
void
NonConformalManager::manage_ghosting(std::vector<stk::mesh::EntityKey>& recvGhostsToRemove)
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  bulk_data.modification_begin();

  if ( nonConformalGhosting_ == NULL) {
    // create new ghosting
    std::string theGhostName = "nalu_nonConformal_ghosting";
    nonConformalGhosting_ = &bulk_data.create_ghosting( theGhostName );
  }
    
  bulk_data.change_ghosting( *nonConformalGhosting_, elemsToGhost_, recvGhostsToRemove);
  
  bulk_data.modification_end();
}

} // namespace nalu
} // namespace sierra
