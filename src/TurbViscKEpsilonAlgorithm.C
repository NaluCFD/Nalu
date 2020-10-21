/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <TurbViscKEpsilonAlgorithm.h>
#include <Algorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// TurbViscKEpsilonAlgorithm - compute tvisc for KEpsilon model
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TurbViscKEpsilonAlgorithm::TurbViscKEpsilonAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part),
    cmu_(realm.get_turb_model_constant(TM_cMu)),
    tke_(NULL),
    density_(NULL),
    tvisc_(NULL)
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  tke_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_ke");
  eps_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_dissipation");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  tvisc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity");
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbViscKEpsilonAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const double small = 1.0e-16;

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*tvisc_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    const double *tke = stk::mesh::field_data(*tke_, *b.begin() );
    const double *eps = stk::mesh::field_data(*eps_, *b.begin() );
    const double *density = stk::mesh::field_data(*density_, *b.begin() );
    double *tvisc = stk::mesh::field_data(*tvisc_, *b.begin() );

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      tvisc[k] = cmu_*density[k]*tke[k]*tke[k]/(eps[k]+small);
    }
  }
}

} // namespace nalu
} // namespace Sierra
