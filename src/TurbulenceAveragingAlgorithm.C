/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <TurbulenceAveragingAlgorithm.h>
#include <Algorithm.h>
#include <AveragingInfo.h>

#include <FieldTypeDef.h>
#include <Realm.h>
#include <TimeIntegrator.h>
#include <NaluEnv.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_io
#include <stk_io/StkMeshIoBroker.hpp>

// c++
#include <vector>
#include <utility>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// TurbulenceAveragingAlgorithm - nodal favre and reynolds averaging
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TurbulenceAveragingAlgorithm::TurbulenceAveragingAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part)
{
  // save off pairs
  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();

  // deal with density; always need Reynolds averaged quantity
  const std::string rhoName = "density";
  const std::string rhoRaName = "density_ra";
  stk::mesh::FieldBase *density = meta_data.get_field(stk::topology::NODE_RANK, rhoName);
  stk::mesh::FieldBase *densityRA = meta_data.get_field(stk::topology::NODE_RANK, rhoRaName);
  reynoldsFieldVecPair_.push_back(std::make_pair(density, densityRA));
  reynoldsFieldSize_.push_back(1);

  // set up Reynolds
  for ( size_t i = 0; i < realm_.averagingInfo_->reynoldsFieldNameVec_.size(); ++i ) {
    std::string primitiveName = realm_.averagingInfo_->reynoldsFieldNameVec_[i];
    std::string averagedName = primitiveName + "_ra";
    stk::mesh::FieldBase *primitive = meta_data.get_field(stk::topology::NODE_RANK, primitiveName);
    stk::mesh::FieldBase *averaged = meta_data.get_field(stk::topology::NODE_RANK, averagedName);
    if ( NULL == primitive || NULL == averaged ) {
      if ( NULL == primitive )
        NaluEnv::self().naluOutputP0() << " Sorry, no primitive field by the name " << primitiveName << std::endl;
      if ( NULL == averaged )
        NaluEnv::self().naluOutputP0() << " Sorry, no reynolds averaged field by the name " << averagedName << std::endl;
      throw std::runtime_error("issue with Reynolds fields existing");
    }
    else {
      reynoldsFieldVecPair_.push_back(std::make_pair(primitive, averaged));
      const unsigned fieldSizePrimitive = primitive->max_size(stk::topology::NODE_RANK);
      const unsigned fieldSizeAveraged = averaged->max_size(stk::topology::NODE_RANK);
      if ( fieldSizePrimitive == fieldSizeAveraged )
        reynoldsFieldSize_.push_back(fieldSizeAveraged);
      else
        throw std::runtime_error("issue with Reynolds field size being equal");
    }
  }

  // set up favre
  for ( size_t i = 0; i < realm_.averagingInfo_->favreFieldNameVec_.size(); ++i ) {
    std::string primitiveName = realm_.averagingInfo_->favreFieldNameVec_[i];
    std::string averagedName = primitiveName + "_fa";
    stk::mesh::FieldBase *primitive = meta_data.get_field(stk::topology::NODE_RANK, primitiveName);
    stk::mesh::FieldBase *averaged = meta_data.get_field(stk::topology::NODE_RANK, averagedName);

    if ( !primitive->type_is<double>() )
      throw std::runtime_error("type is not double");

    if ( NULL == primitive || NULL == averaged) {
      if ( NULL == primitive )
        NaluEnv::self().naluOutputP0() << " Sorry, no primitive field by the name " << primitiveName << std::endl;
      if ( NULL == averaged )
        NaluEnv::self().naluOutputP0() << " Sorry, no favre averaged field by the name " << averagedName << std::endl;
      throw std::runtime_error("issue with Favre fields existing");
    }
    else {
      favreFieldVecPair_.push_back(std::make_pair(primitive, averaged));
      const unsigned fieldSizePrimitive = primitive->max_size(stk::topology::NODE_RANK);
      const unsigned fieldSizeAveraged = averaged->max_size(stk::topology::NODE_RANK);
      if ( fieldSizePrimitive == fieldSizeAveraged )
        favreFieldSize_.push_back(fieldSizeAveraged);
      else
        throw std::runtime_error("issue with Favre field size being equal");
    }
  }

  // review what will be done
  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "Averaging Review:          " << std::endl;
  NaluEnv::self().naluOutputP0() << "=========================" << std::endl;
  for ( size_t iav = 0; iav < reynoldsFieldVecPair_.size(); ++iav ) {
    stk::mesh::FieldBase *primitiveFB = reynoldsFieldVecPair_[iav].first;
    stk::mesh::FieldBase *averageFB = reynoldsFieldVecPair_[iav].second;
    NaluEnv::self().naluOutputP0() << "Primitive/Reynolds name: " << primitiveFB->name() << "/" <<  averageFB->name()
		    << " size " << reynoldsFieldSize_[iav] << std::endl;
  }

  for ( size_t iav = 0; iav < favreFieldVecPair_.size(); ++iav ) {
    stk::mesh::FieldBase *primitiveFB = favreFieldVecPair_[iav].first;
    stk::mesh::FieldBase *averageFB = favreFieldVecPair_[iav].second;
    NaluEnv::self().naluOutputP0() << "Primitive/Favre name:    " << primitiveFB->name() << "/" <<  averageFB->name()
		    << " size " << favreFieldSize_[iav] << std::endl;
  }

}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbulenceAveragingAlgorithm::execute()
{
  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();

  // increment time filter; RESTART for this field...
  const double dt = realm_.timeIntegrator_->get_time_step();
  const double oldTimeFilter = realm_.averagingInfo_->currentTimeFilter_;
  const double timeFilterInterval =  realm_.averagingInfo_->timeFilterInterval_;

  // forced reset
  const bool forcedReset = realm_.averagingInfo_->forcedReset_;

  // check to reset filter
  const bool resetFilter = ( oldTimeFilter + dt  > timeFilterInterval ) || forcedReset;
  const double zeroCurrent = resetFilter ? 0.0 : 1.0;
  const double currentTimeFilter =  resetFilter ? dt : oldTimeFilter + dt;
  realm_.averagingInfo_->currentTimeFilter_ = currentTimeFilter;
  NaluEnv::self().naluOutputP0() << "Filter Size " << currentTimeFilter << std::endl;

  // deactivate hard reset
  realm_.averagingInfo_->forcedReset_ = false;

  // size
  size_t reynoldsFieldPairSize = reynoldsFieldVecPair_.size();
  size_t favreFieldPairSize = favreFieldVecPair_.size();

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // Reynolds averaged density is the first entry
    stk::mesh::FieldBase *densityFB = reynoldsFieldVecPair_[0].first;
    stk::mesh::FieldBase *densityRAFB = reynoldsFieldVecPair_[0].second;
    double *density = (double*)stk::mesh::field_data(*densityFB, b);
    double *densityRA = (double*)stk::mesh::field_data(*densityRAFB, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get node
      stk::mesh::Entity node = b[k];
      
      // save off old density for below Favre procedure
      const double oldRhoRA  = densityRA[k];

      // reynolds first since density is required in Favre
      for ( size_t iav = 0; iav < reynoldsFieldPairSize; ++iav ) {
        stk::mesh::FieldBase *primitiveFB = reynoldsFieldVecPair_[iav].first;
        stk::mesh::FieldBase *averageFB = reynoldsFieldVecPair_[iav].second;
        const double * primitive = (double*)stk::mesh::field_data(*primitiveFB, node);
        double * average = (double*)stk::mesh::field_data(*averageFB, node);
        // get size
        const int fieldSize = reynoldsFieldSize_[iav];
        for ( int j = 0; j < fieldSize; ++j ) {
          const double averageField = (average[j]*oldTimeFilter*zeroCurrent + primitive[j]*dt)/currentTimeFilter;
          average[j] = averageField;
        }
      }

      // save off density for below Favre procedure
      const double rho = density[k];
      const double rhoRA  = densityRA[k];

      // favre
      for ( size_t iav = 0; iav < favreFieldPairSize; ++iav ) {
        stk::mesh::FieldBase *primitiveFB = favreFieldVecPair_[iav].first;
        stk::mesh::FieldBase *averageFB = favreFieldVecPair_[iav].second;
        const double * primitive = (double*)stk::mesh::field_data(*primitiveFB, node);
        double * average = (double*)stk::mesh::field_data(*averageFB,node);
        // get size
        const int fieldSize = favreFieldSize_[iav];
        for ( int j = 0; j < fieldSize; ++j ) {
          const double averageField = (average[j]*oldRhoRA*oldTimeFilter*zeroCurrent + primitive[j]*rho*dt)/currentTimeFilter/rhoRA;
          average[j] = averageField;
        }
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
