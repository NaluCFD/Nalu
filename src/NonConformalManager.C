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
#include <stk_util/environment/CPUTime.hpp>

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
    nonConformalGhosting_(NULL),
    needToGhostCount_(0)
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

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
NonConformalManager::initialize()
{

  const double timeA = stk::cpu_time();

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
 
  // initialize need to ghost and elems to ghost
  needToGhostCount_ = 0;
  elemsToGhost_.clear();

  bulk_data.modification_begin();
  
  if ( nonConformalGhosting_ == NULL) {
    // create new ghosting
    std::string theGhostName = "nalu_nonConformal_ghosting";
    nonConformalGhosting_ = &bulk_data.create_ghosting( theGhostName );
  }
  else {
    bulk_data.destroy_ghosting(*nonConformalGhosting_);
  }
  
  bulk_data.modification_end();
  
  // loop over nonConformalInfo and initialize
  for ( size_t k = 0; k < nonConformalInfoVec_.size(); ++k )
    nonConformalInfoVec_[k]->initialize();
  
  // manage ghosting
  manage_ghosting();
  
  // complete search
  for ( size_t k = 0; k < nonConformalInfoVec_.size(); ++k )
    nonConformalInfoVec_[k]->complete_search();

  // provide diagnosis
  if ( ncAlgDetailedOutput_ ) {
    for ( size_t k = 0; k < nonConformalInfoVec_.size(); ++k )
      nonConformalInfoVec_[k]->provide_diagnosis();
  }

  // end time
  const double timeB = stk::cpu_time();
  realm_.timerContact_ += (timeB-timeA);

}

//--------------------------------------------------------------------------
//-------- manage_ghosting -------------------------------------------------
//--------------------------------------------------------------------------
void
NonConformalManager::manage_ghosting()
{  
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  // check for ghosting need
  uint64_t g_needToGhostCount = 0;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &needToGhostCount_, &g_needToGhostCount, 1);
  if (g_needToGhostCount > 0) {
    
    NaluEnv::self().naluOutputP0() << "NonConformal alg will ghost a number of entities: "
                    << g_needToGhostCount  << std::endl;
    
    bulk_data.modification_begin();
    bulk_data.change_ghosting( *nonConformalGhosting_, elemsToGhost_);
    bulk_data.modification_end();
  }
  else {
    NaluEnv::self().naluOutputP0() << "NonConformal alg will NOT ghost entities: " << std::endl;
  }
}

} // namespace nalu
} // namespace sierra
