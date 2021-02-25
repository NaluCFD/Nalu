/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <SixDofSurfaceForceAndMomentAlgorithmDriver.h>
#include <Algorithm.h>
#include <AlgorithmDriver.h>
#include <FieldFunctions.h>
#include <FieldTypeDef.h>
#include <Realm.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

//==========================================================================
// Class Definition
//==========================================================================
// SurfaceForceAndMomentAlgorithmDriver - Drives six dof integration algs
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SixDofSurfaceForceAndMomentAlgorithmDriver::SixDofSurfaceForceAndMomentAlgorithmDriver( 
  Realm &realm)
  : AlgorithmDriver(realm)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
SixDofSurfaceForceAndMomentAlgorithmDriver::~SixDofSurfaceForceAndMomentAlgorithmDriver()
{
  std::vector<Algorithm *>::iterator iter, iter_end;
  iter_end = algVec_.end();
  for(iter = algVec_.begin(); iter != iter_end; ++iter)
    delete *iter;
}

//--------------------------------------------------------------------------
//-------- zero_fields -------------------------------------------------------
//--------------------------------------------------------------------------
void
SixDofSurfaceForceAndMomentAlgorithmDriver::zero_fields()
{

  // common
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  
  ScalarFieldType *assembledArea = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_area_six_dof");
  // zero fields
  field_fill( meta_data, bulk_data, 0.0, *assembledArea, realm_.get_activate_aura());
}

//--------------------------------------------------------------------------
//-------- parallel_assemble_fields ----------------------------------------
//--------------------------------------------------------------------------
void
SixDofSurfaceForceAndMomentAlgorithmDriver::parallel_assemble_fields()
{

  // Nothing required 

}

//--------------------------------------------------------------------------
//-------- parallel_assemble_area ------------------------------------------
//--------------------------------------------------------------------------
void
SixDofSurfaceForceAndMomentAlgorithmDriver::parallel_assemble_area()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // extract the fields; one of these might be null
  ScalarFieldType *assembledArea = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_area_six_dof");
  // parallel assemble
  std::vector<const stk::mesh::FieldBase*> fields;
  fields.push_back(assembledArea);
  stk::mesh::parallel_sum(bulk_data, fields);

  // periodic assemble
  if ( realm_.hasPeriodic_) {
    const bool bypassFieldCheck = false;
    realm_.periodic_field_update(assembledArea, 1, bypassFieldCheck);
  }

}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
SixDofSurfaceForceAndMomentAlgorithmDriver::execute()
{
  // zero fields
  zero_fields();

  // pre-work; basically, allow algs to load up nodal area(s)
  for ( size_t k = 0; k < algVec_.size(); ++k )
    algVec_[k]->pre_work();

  // assemble area
  parallel_assemble_area();

  // execute
  for ( size_t k = 0; k < algVec_.size(); ++k )
    algVec_[k]->execute();
  
  // parallel assembly
  parallel_assemble_fields();

}


} // namespace nalu
} // namespace Sierra
