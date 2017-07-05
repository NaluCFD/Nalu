/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "overset/UpdateOversetFringeAlgorithmDriver.h"
#include "Realm.h"
#include "overset/OversetManager.h"

#include <stk_mesh/base/Field.hpp>

namespace sierra {
namespace nalu {

UpdateOversetFringeAlgorithmDriver::UpdateOversetFringeAlgorithmDriver(
  Realm& realm)
  : AlgorithmDriver(realm)
{}

UpdateOversetFringeAlgorithmDriver::~UpdateOversetFringeAlgorithmDriver()
{}

void
UpdateOversetFringeAlgorithmDriver::pre_work()
{
  for (auto& f: fields_) {
    realm_.oversetManager_->overset_orphan_node_field_update(
      f->field_, f->sizeRow_, f->sizeCol_);
  }
}

}  // nalu
}  // sierra
