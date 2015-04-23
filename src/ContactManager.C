/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ContactInfo.h>
#include <ContactManager.h>
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
// Acon_ContactManager - manages contact info
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ContactManager::ContactManager(
   Realm &realm)
  : realm_(realm ),
    contactGhosting_(NULL),
    needToGhostCount_(0),
    provideDetailedOutput_(false)
{
  // do nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ContactManager::~ContactManager()
{
  // delete contact info objects:
  std::vector<ContactInfo *>::iterator ic;
  for( ic=contactInfoVec_.begin();
       ic!=contactInfoVec_.end(); ++ic )
    delete (*ic);
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
ContactManager::initialize()
{

  const double timeA = stk::cpu_time();

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
 
  // initialize need to ghost and elems to ghost
  needToGhostCount_ = 0;
  elemsToGhost_.clear();

  bulk_data.modification_begin();
  
  if ( contactGhosting_ == NULL) {
    // create new ghosting
    std::string theGhostName = "nalu_contact_ghosting";
    contactGhosting_ = &bulk_data.create_ghosting( theGhostName );
  }
  else {
    bulk_data.destroy_ghosting(*contactGhosting_);
  }
  
  bulk_data.modification_end();
  
  // loop over contactInfo and initialize
  for ( size_t k = 0; k < contactInfoVec_.size(); ++k )
    contactInfoVec_[k]->initialize();
  
  // manage ghosting
  manage_ghosting();
  
  // complete
  for ( size_t k = 0; k < contactInfoVec_.size(); ++k )
    contactInfoVec_[k]->complete_search();

  // provide output
  if ( provideDetailedOutput_ ) {
    for ( size_t k = 0; k < contactInfoVec_.size(); ++k )
      contactInfoVec_[k]->dump_diagnosis();
  }
  
  // end time
  const double timeB = stk::cpu_time();
  realm_.timerContact_ += (timeB-timeA);

}

//--------------------------------------------------------------------------
//-------- manage_ghosting -------------------------------------------------
//--------------------------------------------------------------------------
void
ContactManager::manage_ghosting()
{
  
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  // check for ghosting need
  uint64_t g_needToGhostCount = 0;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &needToGhostCount_, &g_needToGhostCount, 1);
  if (g_needToGhostCount > 0) {
    
    NaluEnv::self().naluOutputP0() << "Contact alg will ghost a number of entities: "
                    << g_needToGhostCount  << std::endl;
    
    bulk_data.modification_begin();
    bulk_data.change_ghosting( *contactGhosting_, elemsToGhost_);
    bulk_data.modification_end();
  }
  else {
    NaluEnv::self().naluOutputP0() << "Contact alg will NOT ghost entities: " << std::endl;
  }
}


} // namespace nalu
} // namespace sierra
