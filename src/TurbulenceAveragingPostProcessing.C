/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <TurbulenceAveragingPostProcessing.h>
#include <AveragingInfo.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>
#include <Realm.h>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// basic c++
#include <stdexcept>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// TurbulenceAveragingPostProcessing - post process
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TurbulenceAveragingPostProcessing::TurbulenceAveragingPostProcessing(
  Realm & realm,
  const YAML::Node & node) 
  : realm_(realm),
    currentTimeFilter_(0.0),
    timeFilterInterval_(1.0e8),
    forcedReset_(false)
{
  // load the data
  load(node);
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
TurbulenceAveragingPostProcessing::~TurbulenceAveragingPostProcessing()
{
  for ( size_t k = 0; k < averageInfoVec_.size(); ++k )
    delete averageInfoVec_[k];
}

//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbulenceAveragingPostProcessing::load(
  const YAML::Node & y_node)
{
  // output for results
  const YAML::Node *y_average = y_node.FindValue("turbulence_averaging");
  if (y_average) {    
    get_if_present(*y_average, "forced_reset", forcedReset_, forcedReset_);
    get_if_present(*y_average, "time_filter_interval", timeFilterInterval_, timeFilterInterval_);

    // extract the sequence of types
    const YAML::Node *y_specs = expect_sequence(*y_average, "specifications", false);
    if (y_specs) {
      for (size_t ispec = 0; ispec < y_specs->size(); ++ispec) {
        const YAML::Node &y_spec = (*y_specs)[ispec];
        
        // new the info object
        AveragingInfo *avInfo = new AveragingInfo();
        
        // find the name
        const YAML::Node *theName = y_spec.FindValue("name");
        if ( theName )
          *theName >> avInfo->name_;
        else
          throw std::runtime_error("TurbulenceAveragingPostProcessing: no name provided");  
        
        // extract the set of target names
        const YAML::Node &targets = y_spec["target_name"];
        if (targets.Type() == YAML::NodeType::Scalar) {
          avInfo->targetNames_.resize(1);
          targets >> avInfo->targetNames_[0];
        }
        else {
          avInfo->targetNames_.resize(targets.size());
          for (size_t i=0; i < targets.size(); ++i) {
            targets[i] >> avInfo->targetNames_[i];
          }
        }
 
        // reynolds
        const YAML::Node *y_reynolds = y_spec.FindValue("reynolds_averaged_variables");
        if (y_reynolds) {
          size_t varSize = y_reynolds->size();
          for (size_t ioption = 0; ioption < varSize; ++ioption) {
            const YAML::Node & y_var = (*y_reynolds)[ioption];
            std::string fieldName;
            y_var >> fieldName;
            if ( fieldName != "density" )
              avInfo->reynoldsFieldNameVec_.push_back(fieldName);
          }
        }
        
        // favre
        const YAML::Node *y_favre = y_spec.FindValue("favre_averaged_variables");
        if (y_favre) {
          size_t varSize = y_favre->size();
          for (size_t ioption = 0; ioption < varSize; ++ioption) {
            const YAML::Node & y_var = (*y_favre)[ioption];
            std::string fieldName;
            y_var >> fieldName;
            if ( fieldName != "density")
              avInfo->favreFieldNameVec_.push_back(fieldName);
          } 
        }
        
        // check for reynolds stress and tke post processing
        get_if_present(y_spec, "compute_reynolds_stress", avInfo->computeReynoldsStress_, avInfo->computeReynoldsStress_);
        get_if_present(y_spec, "compute_tke", avInfo->computeTke_, avInfo->computeTke_);
        
        // we will need Reynolds-averaged velocity if we need to compute TKE
        if ( avInfo->computeTke_ || avInfo->computeReynoldsStress_ ) {
          const std::string velocityName = "velocity";
          if ( std::find(avInfo->reynoldsFieldNameVec_.begin(), avInfo->reynoldsFieldNameVec_.end(), velocityName) == avInfo->reynoldsFieldNameVec_.end() ) {
            // not found; add it
            avInfo->reynoldsFieldNameVec_.push_back(velocityName);
          }  
        } 
        
        // push back the object
        averageInfoVec_.push_back(avInfo);
      }
    }
    else {
      throw std::runtime_error("TurbulenceAveragingPostProcessing: no specifications provided");
    }
  }
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbulenceAveragingPostProcessing::setup()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();

  // loop over all info and setup (register fields, set parts, etc.)
  for (size_t k = 0; k < averageInfoVec_.size(); ++k ) {
 
    // extract the turb info and the name
    AveragingInfo *avInfo = averageInfoVec_[k];

    const std::string averageBlockName = avInfo->name_;

    // loop over all target names, extract the part; register the fields
    for ( size_t itarget = 0; itarget < avInfo->targetNames_.size(); ++itarget ) {
      stk::mesh::Part *targetPart = metaData.get_part(avInfo->targetNames_[itarget]);
      if ( NULL == targetPart ) {
        NaluEnv::self().naluOutputP0() << "Trouble with part " << avInfo->targetNames_[itarget] << std::endl;
        throw std::runtime_error("Sorry, no part name found by the name: " + avInfo->targetNames_[itarget]);
      }
      else {
        // push back
        avInfo->partVec_.push_back(targetPart);
      }
      
      // first, register TKE
      if ( avInfo->computeTke_ ) {
        // hack a name; the name is not tied to the average info name
        const std::string tkeName = "resolved_turbulent_ke";
        // register and put the field
        stk::mesh::FieldBase *tkeField 
          = &(metaData.declare_field< stk::mesh::Field<double, stk::mesh::SimpleArrayTag> >(stk::topology::NODE_RANK, tkeName));
        stk::mesh::put_field(*tkeField,*targetPart,1);
        // augment the restart list
        realm_.augment_restart_variable_list(tkeName);
      }

      // second, register stress
      if ( avInfo->computeReynoldsStress_ ) {
        // hack a name; the name is not tied to the average info name
        const std::string stressName = "reynolds_stress";
        // register and put the field
        stk::mesh::FieldBase *stressField 
          = &(metaData.declare_field< stk::mesh::Field<double, stk::mesh::SimpleArrayTag> >(stk::topology::NODE_RANK, stressName));
        // only output the unique components of the tensor
        const int stressSize = realm_.spatialDimension_ == 3 ? 6 : 3;
        stk::mesh::put_field(*stressField, *targetPart, stressSize);
        // augment the restart list
        realm_.augment_restart_variable_list(stressName);
      }

      // deal with density; always need Reynolds averaged quantity
      const std::string densityReynoldsName = "density_ra_" + averageBlockName;
      ScalarFieldType *densityReynolds =  &(metaData.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, densityReynoldsName));
      stk::mesh::put_field(*densityReynolds, *targetPart);
      
      // reynolds
      for ( size_t i = 0; i < avInfo->reynoldsFieldNameVec_.size(); ++i ) {
        const std::string primitiveName = avInfo->reynoldsFieldNameVec_[i];
        const std::string averagedName = primitiveName + "_ra_" + averageBlockName;
        register_field(primitiveName, averagedName, metaData, targetPart);
      }
      
      // favre
      for ( size_t i = 0; i < avInfo->favreFieldNameVec_.size(); ++i ) {
        const std::string primitiveName = avInfo->favreFieldNameVec_[i];
        const std::string averagedName = primitiveName + "_fa_" + averageBlockName;
        register_field(primitiveName, averagedName, metaData, targetPart);
      }      
    }

    // now deal with pairs; extract density
    const std::string densityName = "density";
    const std::string densityReynoldsName = "density_ra_" + averageBlockName;
    stk::mesh::FieldBase *density = metaData.get_field(stk::topology::NODE_RANK, densityName);
    stk::mesh::FieldBase *densityReynolds = metaData.get_field(stk::topology::NODE_RANK, densityReynoldsName);
    avInfo->reynoldsFieldVecPair_.push_back(std::make_pair(density, densityReynolds));
    avInfo->reynoldsFieldSizeVec_.push_back(1);
    realm_.augment_restart_variable_list(densityReynoldsName);
    
    // reynolds
    for ( size_t i = 0; i < avInfo->reynoldsFieldNameVec_.size(); ++i ) {
      const std::string primitiveName = avInfo->reynoldsFieldNameVec_[i];
      const std::string averagedName = primitiveName + "_ra_" + averageBlockName;
      construct_pair(primitiveName, averagedName, avInfo->reynoldsFieldVecPair_, avInfo->reynoldsFieldSizeVec_, metaData);
    }
    
    // favre
    for ( size_t i = 0; i < avInfo->favreFieldNameVec_.size(); ++i ) {
      const std::string primitiveName = avInfo->favreFieldNameVec_[i];
      const std::string averagedName = primitiveName + "_fa_" + averageBlockName;
      construct_pair(primitiveName, averagedName, avInfo->favreFieldVecPair_, avInfo->favreFieldSizeVec_, metaData);
    }

    // output what we have done here...
    review(avInfo);
  }
}

//--------------------------------------------------------------------------
//-------- register_field --------------------------------------------------
//--------------------------------------------------------------------------
void
TurbulenceAveragingPostProcessing::register_field(
  const std::string primitiveName,
  const std::string averagedName,
  stk::mesh::MetaData &metaData,
  stk::mesh::Part *part)
{
  // first, augment the restart list
  realm_.augment_restart_variable_list(averagedName);

  // declare field; put the field and augment restart; need size from the primitive
  stk::mesh::FieldBase *primitiveField = metaData.get_field(stk::topology::NODE_RANK, primitiveName);

  // check for existance and if it is a double
  if ( NULL == primitiveField )
    throw std::runtime_error("TurbulenceAveragingPostProcessing::register_field() no primitive by this name: " + primitiveName);

  if ( !primitiveField->type_is<double>() )
    throw std::runtime_error("TurbulenceAveragingPostProcessing::register_field() type of field is not double: " +primitiveName);
  
  // extract size (would love to do this by part), however, not yet a use case
  const unsigned fieldSizePrimitive = primitiveField->max_size(stk::topology::NODE_RANK);

  // register the averaged field with this size; treat velocity as a special case to retain the vector aspect
  if ( primitiveName == "velocity" ) {
    VectorFieldType *averagedField = &(metaData.declare_field<VectorFieldType>(stk::topology::NODE_RANK, averagedName));
    stk::mesh::put_field(*averagedField, *part, fieldSizePrimitive);
  }
  else {
    stk::mesh::FieldBase *averagedField 
      = &(metaData.declare_field< stk::mesh::Field<double, stk::mesh::SimpleArrayTag> >(stk::topology::NODE_RANK, averagedName));
    stk::mesh::put_field(*averagedField, *part, fieldSizePrimitive);
  }
}
  
//--------------------------------------------------------------------------
//-------- construct_pair --------------------------------------------------
//--------------------------------------------------------------------------
void
TurbulenceAveragingPostProcessing::construct_pair(
  const std::string primitiveName,
  const std::string averagedName,
  std::vector<std::pair<stk::mesh::FieldBase *, stk::mesh::FieldBase *> > &fieldVecPair,
  std::vector<unsigned> &fieldSizeVec,
  stk::mesh::MetaData &metaData)
{ 
  // augment the restart list
  realm_.augment_restart_variable_list(averagedName);

  // extract the valid primitive and averaged field
  stk::mesh::FieldBase *primitiveField = metaData.get_field(stk::topology::NODE_RANK, primitiveName);
  stk::mesh::FieldBase *averagedField = metaData.get_field(stk::topology::NODE_RANK, averagedName);
  
  // the size; gaurenteed to be the same based on the field registration
  const unsigned fieldSizeAveraged = averagedField->max_size(stk::topology::NODE_RANK);
  fieldSizeVec.push_back(fieldSizeAveraged);
  
  // construct pairs
  fieldVecPair.push_back(std::make_pair(primitiveField, averagedField));
}

//--------------------------------------------------------------------------
//-------- review ----------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbulenceAveragingPostProcessing::review( 
  const AveragingInfo *avInfo)
{
  // review what will be done
  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "Averaging Review: " << avInfo->name_ << std::endl;
  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
  for ( size_t iav = 0; iav < avInfo->reynoldsFieldVecPair_.size(); ++iav ) {
    stk::mesh::FieldBase *primitiveFB = avInfo->reynoldsFieldVecPair_[iav].first;
    stk::mesh::FieldBase *averageFB = avInfo->reynoldsFieldVecPair_[iav].second;
    NaluEnv::self().naluOutputP0() << "Primitive/Reynolds name: " << primitiveFB->name() << "/" <<  averageFB->name()
                                   << " size " << avInfo->reynoldsFieldSizeVec_[iav] << std::endl;
  }
  
  for ( size_t iav = 0; iav < avInfo->favreFieldVecPair_.size(); ++iav ) {
    stk::mesh::FieldBase *primitiveFB = avInfo->favreFieldVecPair_[iav].first;
    stk::mesh::FieldBase *averageFB = avInfo->favreFieldVecPair_[iav].second;
    NaluEnv::self().naluOutputP0() << "Primitive/Favre name:    " << primitiveFB->name() << "/" <<  averageFB->name()
                                   << " size " << avInfo->favreFieldSizeVec_[iav] << std::endl;
  }
  
  if ( avInfo->computeTke_ ) {
    NaluEnv::self().naluOutputP0() << "TKE will be computed; add resolved_turbulent_ke to the Reynolds/Favre block for mean"<< std::endl;
  }
  
  if ( avInfo->computeReynoldsStress_ ) {
    NaluEnv::self().naluOutputP0() << "Reynolds Stress will be computed; add reynolds_stress to output"<< std::endl;
  }
  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbulenceAveragingPostProcessing::execute()
{
  stk::mesh::MetaData &metaData = realm_.meta_data();

  // increment time filter; RESTART for this field...
  const double dt = realm_.get_time_step();
  const double oldTimeFilter = currentTimeFilter_;

  // check to reset filter
  const bool resetFilter = ( oldTimeFilter + dt  > timeFilterInterval_ ) || forcedReset_;
  const double zeroCurrent = resetFilter ? 0.0 : 1.0;
  currentTimeFilter_ =  resetFilter ? dt : oldTimeFilter + dt;
  NaluEnv::self().naluOutputP0() << "Filter Size " << currentTimeFilter_ << std::endl;

  // deactivate hard reset
  forcedReset_ = false;

  // loop over all info and setup (register fields, set parts, etc.)
  for (size_t k = 0; k < averageInfoVec_.size(); ++k ) {

    // extract the turb info and the name
    AveragingInfo *avInfo = averageInfoVec_[k];

    // size
    size_t reynoldsFieldPairSize = avInfo->reynoldsFieldVecPair_.size();
    size_t favreFieldPairSize = avInfo->favreFieldVecPair_.size();

    // define some common selectors
    stk::mesh::Selector s_all_nodes
      = (metaData.locally_owned_part() | metaData.globally_shared_part())
      & stk::mesh::selectUnion(avInfo->partVec_) 
      & !(realm_.get_inactive_selector());

    stk::mesh::BucketVector const& node_buckets =
      realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
    for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
          ib != node_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();
      
      // Reynolds averaged density is the first entry (FieldBase == FB)
      stk::mesh::FieldBase *densityFB = avInfo->reynoldsFieldVecPair_[0].first;
      stk::mesh::FieldBase *densityRAFB = avInfo->reynoldsFieldVecPair_[0].second;
      double *density = (double*)stk::mesh::field_data(*densityFB, b);
      double *densityRA = (double*)stk::mesh::field_data(*densityRAFB, b);
      
      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

        // get node
        stk::mesh::Entity node = b[k];
      
        // save off old density for below Favre procedure
        const double oldRhoRA  = densityRA[k];
        
        // reynolds first since density is required in Favre
        for ( size_t iav = 0; iav < reynoldsFieldPairSize; ++iav ) {
          stk::mesh::FieldBase *primitiveFB = avInfo->reynoldsFieldVecPair_[iav].first;
          stk::mesh::FieldBase *averageFB = avInfo->reynoldsFieldVecPair_[iav].second;
          const double * primitive = (double*)stk::mesh::field_data(*primitiveFB, node);
          double * average = (double*)stk::mesh::field_data(*averageFB, node);
          // get size
          const int fieldSize = avInfo->reynoldsFieldSizeVec_[iav];
          for ( int j = 0; j < fieldSize; ++j ) {
            const double averageField = (average[j]*oldTimeFilter*zeroCurrent + primitive[j]*dt)/currentTimeFilter_;
            average[j] = averageField;
          }
        }

        // save off density for below Favre procedure
        const double rho = density[k];
        const double rhoRA  = densityRA[k];

        // favre
        for ( size_t iav = 0; iav < favreFieldPairSize; ++iav ) {
          stk::mesh::FieldBase *primitiveFB = avInfo->favreFieldVecPair_[iav].first;
          stk::mesh::FieldBase *averageFB = avInfo->favreFieldVecPair_[iav].second;
          const double * primitive = (double*)stk::mesh::field_data(*primitiveFB, node);
          double * average = (double*)stk::mesh::field_data(*averageFB,node);
          // get size
          const int fieldSize = avInfo->favreFieldSizeVec_[iav];
          for ( int j = 0; j < fieldSize; ++j ) {
            const double averageField = (average[j]*oldRhoRA*oldTimeFilter*zeroCurrent + primitive[j]*rho*dt)/currentTimeFilter_/rhoRA;
            average[j] = averageField;
          }
        }
      }
    }
  
    // process tke
    if ( avInfo->computeTke_ ) {

      const int nDim = realm_.spatialDimension_;

      const std::string averageBlockName = avInfo->name_;
      const std::string velocityRAName = "velocity_ra_" + averageBlockName;
      const std::string resolvedTkeName = "resolved_turbulent_ke";

      // extract fields
      stk::mesh::FieldBase *velocity = metaData.get_field(stk::topology::NODE_RANK, "velocity");
      stk::mesh::FieldBase *velocityRA = metaData.get_field(stk::topology::NODE_RANK, velocityRAName);
      stk::mesh::FieldBase *resololvedTke = metaData.get_field(stk::topology::NODE_RANK, resolvedTkeName);
      
      stk::mesh::BucketVector const& node_buckets_tke =
        realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
      for ( stk::mesh::BucketVector::const_iterator ib = node_buckets_tke.begin();
            ib != node_buckets_tke.end() ; ++ib ) {
        stk::mesh::Bucket & b = **ib ;
        const stk::mesh::Bucket::size_type length   = b.size();
      
        // fields
        const double *uNp1 = (double*)stk::mesh::field_data(*velocity, b);
        const double *uNp1RA = (double*)stk::mesh::field_data(*velocityRA, b);
        double *tke = (double*)stk::mesh::field_data(*resololvedTke, b);
        
        for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
          double sum = 0.0;
          for ( int j = 0; j < nDim; ++j ) {
            const double uPrime = uNp1[k*nDim+j] - uNp1RA[k*nDim+j];
            sum += 0.5*uPrime*uPrime;
          }
          tke[k] = sum;
        }
      }
    } 
  
    // process stress
    if ( avInfo->computeReynoldsStress_ ) {

      const int nDim = realm_.spatialDimension_;
      const int stressSize = realm_.spatialDimension_ == 3 ? 6 : 3;
      
      const std::string averageBlockName = avInfo->name_;
      const std::string velocityRAName = "velocity_ra_" + averageBlockName;
      const std::string stressName = "reynolds_stress";

      // extract fields
      stk::mesh::FieldBase *velocity = metaData.get_field(stk::topology::NODE_RANK, "velocity");
      stk::mesh::FieldBase *velocityRA = metaData.get_field(stk::topology::NODE_RANK, velocityRAName);
      stk::mesh::FieldBase *reynoldsStress = metaData.get_field(stk::topology::NODE_RANK, stressName);

      stk::mesh::BucketVector const& node_buckets_stress =
        realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
      for ( stk::mesh::BucketVector::const_iterator ib = node_buckets_stress.begin();
            ib != node_buckets_stress.end() ; ++ib ) {
        stk::mesh::Bucket & b = **ib ;
        const stk::mesh::Bucket::size_type length   = b.size();   

        // fields
        const double *uNp1 = (double*)stk::mesh::field_data(*velocity, b);
        const double *uNp1RA = (double*)stk::mesh::field_data(*velocityRA, b);
        double *stress = (double*)stk::mesh::field_data(*reynoldsStress, b);
        
        for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
          
          // stress is symmetric, so only save off 6 or 3 components
          int componentCount = 0;
          for ( int i = 0; i < nDim; ++i ) {
            const double ui = uNp1[k*nDim+i];
            const double uiRA = uNp1RA[k*nDim+i];
            for ( int j = i; j < nDim; ++j ) {
              const int component = componentCount;
              const double uj = uNp1[k*nDim+j];
              const double ujRA = uNp1RA[k*nDim+j];
              const double newStress = (stress[k*stressSize+component]*oldTimeFilter*zeroCurrent + ui*uj*dt - uiRA*ujRA*dt)/currentTimeFilter_;
              stress[k*stressSize+component] = newStress;
              componentCount++;
            }
          }
        }
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
