/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <BurkeSchumannAlgorithm.h>
#include <Algorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <ReferencePropertyData.h>
#include <stk_mesh/base/Field.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_io
#include <stk_io/StkMeshIoBroker.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// BurkeSchumannAlgorithm - post process species from Z
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
BurkeSchumannAlgorithm::BurkeSchumannAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap)
  : Algorithm(realm, part),
    massFracSize_(5),
    mixtureFraction_(NULL),
    massFraction_(NULL),
    zStoich_(0.0),
    fuelId_(0),
    oxidizerId_(1),
    carbonDioxideId_(2),
    waterId_(3),
    diluentId_(4),
    M_(1.0), /* hard coded for CH4 and O2 combustion */
    N_(4.0)
{

  // save off fields
  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();
  mixtureFraction_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixture_fraction");
  massFraction_ = meta_data.get_field<GenericFieldType>(stk::topology::NODE_RANK, "mass_fraction");

  // size vectors
  mwVec_.resize(massFracSize_);
  stoichiometryVec_.resize(massFracSize_);
  primaryVec_.resize(massFracSize_);
  secondaryVec_.resize(massFracSize_);

  // fill internal vectors
  size_t k = 0;
  std::map<std::string, ReferencePropertyData*>::const_iterator itrp;
  for ( itrp = referencePropertyDataMap.begin();
        itrp!= referencePropertyDataMap.end(); ++itrp, ++k) {
    ReferencePropertyData *propData = (*itrp).second;
    mwVec_[k] = propData->mw_;
    stoichiometryVec_[k] = propData->stoichiometry_;
    primaryVec_[k] = propData->primaryMassFraction_;
    secondaryVec_[k] = propData->secondaryMassFraction_;
  }

  // compute stoichiometric mixture fraction
  zStoich_ = 1.0/(1.0+(M_+N_/4.0)*primaryVec_[fuelId_]*mwVec_[oxidizerId_]/(secondaryVec_[oxidizerId_]*mwVec_[fuelId_]));
  Env::outputP0() << "Stoichiometric Mixture Fraction " << zStoich_ << std::endl;
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
BurkeSchumannAlgorithm::execute()
{

  const bool doCombustion = true;

  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();

  // deal with state
  ScalarFieldType &mixFracNp1 = mixtureFraction_->field_of_state(stk::mesh::StateNP1);

  // pointers and size
  const int nm1Size = massFracSize_-1;
  const double *pt_mw = &mwVec_[0];
  const double *pt_primaryVec = &primaryVec_[0];
  const double *pt_secondaryVec = &secondaryVec_[0];

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*mixtureFraction_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );

  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
	ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    const double *mixFrac = stk::mesh::field_data(mixFracNp1, b);
    double *massFrac = stk::mesh::field_data(*massFraction_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      const double z = mixFrac[k];

      // stoichiometric composition for products
      if ( doCombustion ) {
        const double yCarbonDioxideStoich = pt_primaryVec[fuelId_]*zStoich_*M_*pt_mw[carbonDioxideId_]/pt_mw[fuelId_];
        const double yWaterStoich = pt_primaryVec[fuelId_]*zStoich_*N_*pt_mw[waterId_]/pt_mw[fuelId_]/2.0;
        if ( z <= zStoich_ ) {
          // fuel all consumed
          massFrac[k*massFracSize_ + fuelId_] = 0.0;
          massFrac[k*massFracSize_ + oxidizerId_] = pt_secondaryVec[oxidizerId_]*(1.0-z/zStoich_);
          // products
          const double fac = z/zStoich_;
          massFrac[k*massFracSize_ + carbonDioxideId_] = yCarbonDioxideStoich*fac;
          massFrac[k*massFracSize_ + waterId_] = yWaterStoich*fac;
        }
        else {
          // oxygen all consumed
          massFrac[k*massFracSize_ + oxidizerId_] = 0.0;
          massFrac[k*massFracSize_ + fuelId_] = pt_primaryVec[fuelId_]*(z-zStoich_)/(1.0-zStoich_);
          // products
          const double fac = (1.0-z)/(1.0-zStoich_);
          massFrac[k*massFracSize_ + carbonDioxideId_] = yCarbonDioxideStoich*fac;
          massFrac[k*massFracSize_ + waterId_] = yWaterStoich*fac;
        }

        // n-1 diluent
        double sum = 0.0;
        for ( int j = 0; j < nm1Size; ++j ) {
          sum += massFrac[k*massFracSize_ +j];
        }
        massFrac[k*massFracSize_ + diluentId_] = 1.0-sum;
      }
      else {
	for ( size_t j = 0; j < massFracSize_; ++j ) {
          massFrac[k*massFracSize_ +j] = z*pt_primaryVec[j] + (1.0-z)*pt_secondaryVec[j];
	}
	/*double nm1Sum = 0.0;
        for ( int j = 0; j < nm1Size; ++j ) {
	  const double yj = z*pt_primaryVec[j] + (1.0-z)*pt_secondaryVec[j];
	  const double yjClipped = std::min(1.0, std::max(0.0, yj));
          massFrac[k*massFracSize_ +j] = yjClipped;
	  nm1Sum += yjClipped;
        }
	massFrac[k*massFracSize_+diluentId_] = 1.0-nm1Sum;
	*/
      }

    }

  }
}

} // namespace nalu
} // namespace Sierra
