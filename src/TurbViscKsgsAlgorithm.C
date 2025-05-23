/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <TurbViscKsgsAlgorithm.h>
#include <Algorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SolutionOptions.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Field.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// TurbViscKsgsAlgorithm - compute tvisc for Ksgs model
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TurbViscKsgsAlgorithm::TurbViscKsgsAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part),
    tke_(NULL),
    density_(NULL),
    tvisc_(NULL),
    dualNodalVolume_(NULL),
    cmuEps_(NULL),
    cEps_(NULL),
    visc_(NULL),
    minDistance_(NULL),
    dsqrtkSq_(NULL),
    lrksgsfac_(0.0),
    Cl_(4.0),
    Ao_(20.0),
    Bo_(2.0)
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  tke_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "turbulent_ke");
  density_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "density");
  tvisc_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "turbulent_viscosity");
  dualNodalVolume_ =meta_data.get_field<double>(stk::topology::NODE_RANK, "dual_nodal_volume");
  cmuEps_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "c_mu_epsilon");

  // low-Re form
  cEps_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "c_epsilon");
  visc_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "viscosity");
  // assign required variables that may not be registered to an arbitrary field
  minDistance_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "viscosity");
  dsqrtkSq_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "viscosity");
  if (realm_.solutionOptions_->turbulenceModel_ == LRKSGS ) {
    minDistance_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "minimum_distance_to_wall");
    dsqrtkSq_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "dsqrtk_dx_sq");
    lrksgsfac_ = 1.0;
  }
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbViscKsgsAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const double invNdim = 1.0/meta_data.spatial_dimension();

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
    const double *density = stk::mesh::field_data(*density_, *b.begin() );
    const double *dualNodalVolume = stk::mesh::field_data(*dualNodalVolume_, *b.begin() );
    const double *cmuEps = stk::mesh::field_data(*cmuEps_, *b.begin() );
    // low-Re
    const double *cEps = stk::mesh::field_data(*cEps_, *b.begin() );
    const double *visc = stk::mesh::field_data(*visc_, *b.begin() );
    const double *minDistance = stk::mesh::field_data(*minDistance_, *b.begin() );
    const double *dsqrtkSq = stk::mesh::field_data(*dsqrtkSq_, *b.begin() );
    double *tvisc = stk::mesh::field_data(*tvisc_, *b.begin() );

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const double filter = std::pow(dualNodalVolume[k], invNdim);
      // low-Re corrections
      const double dsq = lrksgsfac_*dsqrtkSq[k];
      const double md = lrksgsfac_*minDistance[k];
      const double nu = visc[k]/density[k];
      const double epsKsgs = lrksgsfac_*(cEps[k]*std::pow(tke[k], 1.5)/filter + 2.0*nu*dsq);
      const double upeps = std::pow(nu*epsKsgs, 0.25)*std::sqrt(Cl_*md/filter);
      const double ypeps = md*upeps/nu;
      const double fmu = lrksgsfac_*(1.0 - std::exp(-std::pow(std::pow(ypeps/Ao_, 2.0/3.0),Bo_))) 
        + (1.0-lrksgsfac_);
      tvisc[k] = fmu*cmuEps[k]*density[k]*std::sqrt(tke[k])*filter;
    }
  }
}
  
} // namespace nalu
} // namespace Sierra
