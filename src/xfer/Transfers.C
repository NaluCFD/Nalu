/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <xfer/Transfers.h>
#include <xfer/Transfer.h>
#include <Simulation.h>
#include <Realms.h>
#include <Realm.h>

// yaml for parsing..
#include <yaml-cpp/yaml.h>

#include <stk_mesh/base/BulkData.hpp>

// basic c++
#include <vector>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// Transfers - do some stuff
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Transfers::Transfers( 
  Simulation &sim)
  : simulation_(sim)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Transfers::~Transfers()
{
  for (size_t ir = 0; ir < transferVector_.size(); ir++)
    delete transferVector_[ir];
}

void 
Transfers::load(const YAML::Node & node) 
{
  // xfers are optional...
  const YAML::Node *transfers = node.FindValue("transfers");
  if (transfers) {
    for ( size_t itransfer = 0; itransfer < transfers->size(); ++itransfer ) {
      const YAML::Node & transferNode = (*transfers)[itransfer];
      Transfer *transferInfo = new Transfer(*this);
      transferInfo->load(transferNode);
      transferVector_.push_back(transferInfo);
    }
  }
}
  
void 
Transfers::breadboard()
{
  for ( size_t itransfer = 0; itransfer < transferVector_.size(); ++itransfer ) {
    transferVector_[itransfer]->breadboard();
  }
}

void 
Transfers::initialize()
{
  for ( size_t itransfer = 0; itransfer < transferVector_.size(); ++itransfer ) {
    transferVector_[itransfer]->initialize_begin();
  }

  for ( size_t itransfer = 0; itransfer < transferVector_.size(); ++itransfer ) {
    const std::string fromName = transferVector_[itransfer]->realmPairName_.first;
    stk::mesh::BulkData &fromBulkData = root()->realms_->find_realm(fromName)->bulk_data();
    fromBulkData.modification_begin();
    transferVector_[itransfer]->change_ghosting(); 
    fromBulkData.modification_end();
  }

  for ( size_t itransfer = 0; itransfer < transferVector_.size(); ++itransfer ) {
    transferVector_[itransfer]->initialize_end();
  }
}

Simulation *Transfers::root() { return parent()->root(); }
Simulation *Transfers::parent() { return &simulation_; }

} // namespace nalu
} // namespace Sierra
