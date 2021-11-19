/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <gas_dynamics/AssembleGasDynamicsAlgorithmDriver.h>
#include <Algorithm.h>
#include <AlgorithmDriver.h>
#include <FieldTypeDef.h>
#include <FieldFunctions.h>
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
// AssembleGasDynamicsAlgorithmDriver - Drives RHS assembly
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleGasDynamicsAlgorithmDriver::AssembleGasDynamicsAlgorithmDriver(
  Realm &realm)
  : AlgorithmDriver(realm)
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleGasDynamicsAlgorithmDriver::~AssembleGasDynamicsAlgorithmDriver()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- pre_work --------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleGasDynamicsAlgorithmDriver::pre_work()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  // zero
  GenericFieldType *rhsGasDyn = metaData.get_field<GenericFieldType>(stk::topology::NODE_RANK, "rhs_gas_dynamics");
  field_fill( metaData, bulkData, 0.0, *rhsGasDyn, realm_.get_activate_aura());
}

//--------------------------------------------------------------------------
//-------- post_work -------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleGasDynamicsAlgorithmDriver::post_work()
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  stk::mesh::MetaData & metaData = realm_.meta_data();

  GenericFieldType *rhsGasDyn = metaData.get_field<GenericFieldType>(stk::topology::NODE_RANK, "rhs_gas_dynamics");

  // u, v, w + cont + e 
  const unsigned totalSize  = metaData.spatial_dimension() + 2;

  // parallel sum
  std::vector<const stk::mesh::FieldBase*> sum_fields(1, rhsGasDyn);
  stk::mesh::parallel_sum(bulkData, sum_fields);

  // periodic
  if ( realm_.hasPeriodic_) {
    realm_.periodic_field_update(rhsGasDyn, totalSize);
  }

  // overset
  if ( realm_.hasOverset_ ) {
    realm_.overset_constraint_node_field_update(rhsGasDyn, 1, totalSize);
  }
}

} // namespace nalu
} // namespace Sierra
