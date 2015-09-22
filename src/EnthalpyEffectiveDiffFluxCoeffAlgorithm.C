/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <EnthalpyEffectiveDiffFluxCoeffAlgorithm.h>
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
// EnthalpyEffectiveDiffFluxCoeffAlgorithm - compute h eff. diff flux coeff
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyEffectiveDiffFluxCoeffAlgorithm::EnthalpyEffectiveDiffFluxCoeffAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *thermalCond,
  ScalarFieldType *specHeat,
  ScalarFieldType *tvisc,
  ScalarFieldType *evisc,
  const double sigmaTurb)
  : Algorithm(realm, part),
    thermalCond_(thermalCond),
    specHeat_(specHeat),
    tvisc_(tvisc),
    evisc_(evisc),
    sigmaTurb_(sigmaTurb),
    isTurbulent_(realm_.is_turbulent())
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
EnthalpyEffectiveDiffFluxCoeffAlgorithm::execute()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*specHeat_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );

  if ( isTurbulent_ ) {
    const double invSigmaTurb = 1.0/sigmaTurb_;

    for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
          ib != node_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();

      const double * thermalCond = stk::mesh::field_data(*thermalCond_, b);
      const double * specHeat = stk::mesh::field_data(*specHeat_, b);
      const double * tvisc = stk::mesh::field_data(*tvisc_, b);
      double * evisc = stk::mesh::field_data(*evisc_, b);

      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        evisc[k] = thermalCond[k]/specHeat[k] + tvisc[k]*invSigmaTurb;
      }
    }
  }
  else {
    for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
          ib != node_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();

      const double * thermalCond = stk::mesh::field_data(*thermalCond_, b);
      const double * specHeat = stk::mesh::field_data(*specHeat_, b);
      double * evisc = stk::mesh::field_data(*evisc_, b);

      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        evisc[k] = thermalCond[k]/specHeat[k];
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
