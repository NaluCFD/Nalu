/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ExplicitFiltering.h>
#include <FieldFunctions.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>
#include <NaluEnv.h>
#include <Realm.h>
#include <SolutionOptions.h>

// master elements
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

// stk_search
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>

// basic c++
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <cmath>

namespace sierra{
namespace nalu{


//==========================================================================
// Class Definition
//==========================================================================
// ExplicitFiltering - assemble a box filter to a node
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ExplicitFiltering::ExplicitFiltering(
  Realm &realm,
  const YAML::Node &node)
  : realm_(realm),
    searchMethod_(stk::search::KDTREE),
    explicitFilteringGhosting_(NULL),
    needToGhostCount_(0),
    debugOutput_(false),
    normalizeResidual_(false)
{
  // load the data
  load(node);
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ExplicitFiltering::~ExplicitFiltering()
{
  // nothing to delete
}


//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
ExplicitFiltering::load(
  const YAML::Node & y_node)
{
  // check for any data probes
  const YAML::Node y_explicit = y_node["explicit_filtering"];
  if (y_explicit) {
    NaluEnv::self().naluOutputP0() << "ExplicitFiltering::load" << std::endl;
    
    // extract the set of from target names; each spec is homogeneous in this respect
    const YAML::Node searchTargets = y_explicit["search_target_part"];
    if (searchTargets.Type() == YAML::NodeType::Scalar) {
      searchTargetNames_.resize(1);
      searchTargetNames_[0] = searchTargets.as<std::string>() ;
    }
    else {
      searchTargetNames_.resize(searchTargets.size());
      for (size_t i=0; i < searchTargets.size(); ++i) {
        searchTargetNames_[i] = searchTargets[i].as<std::string>() ;
      }
    }

    const YAML::Node filterSize = y_explicit["filter_size"];
    if ( filterSize )
      filterSize_ = filterSize.as<Coordinates>() ;
    else
      throw std::runtime_error("ExplicitFiltering::error: filter_size must be provided");
    
    const YAML::Node debugOutput = y_explicit["debug_output"];
    if ( debugOutput ) {
      debugOutput_ = debugOutput.as<bool>();
      NaluEnv::self().naluOutputP0() << "ExplicitFiltering::debugOutput_: ACTIVE" << std::endl;
    }
    
    // extract information from this specification
    const YAML::Node normalizeResidual = y_explicit["normalize_residual"];
    if ( normalizeResidual ) {
      normalizeResidual_ = normalizeResidual.as<bool>();
      NaluEnv::self().naluOutputP0() << "Normalization active? "  << normalizeResidual_ << std::endl;
    }

  }

  // extract the sequence of types
  const YAML::Node y_especs = expect_sequence(y_explicit, "explicit_specifications", true);
  if (y_especs) {
    for (size_t ispec = 0; ispec < y_especs.size(); ++ispec) {
      const YAML::Node y_espec = y_especs[ispec] ;
      
      // find the name, size and type
      const YAML::Node fieldNameNode = y_espec["field_name"];
      const YAML::Node expFieldNameNode = y_espec["explicit_field_name"];
      const YAML::Node fieldSizeNode = y_espec["field_size"];
      const YAML::Node fieldTypeNode = y_espec["field_type"];
      const YAML::Node fieldStateSizeNode = y_espec["field_state_size"];
      
      if ( !fieldNameNode ) 
        throw std::runtime_error("ExplicitFiltering::load():error(), field_name must be provided");
      
      if ( !expFieldNameNode ) 
        throw std::runtime_error("ExplicitFiltering::load():error(), explicit_field_name must be provided");
      
      if ( !fieldSizeNode ) 
        throw std::runtime_error("ExplicitFiltering::load():error(), field_size must be provided");
      
      if ( fieldTypeNode ) {
        std::string fieldType = fieldTypeNode.as<std::string>();
        if ( fieldType != "node_rank" )
          throw std::runtime_error("ExplicitFiltering::load():error() only supports nodal_rank types");
      }

      // extract field state size; optional with a default of unity
      int fieldStateSize = 1;
      if ( fieldStateSizeNode ) {
        fieldStateSize = fieldStateSizeNode.as<int>();
      }
      
      // create the struct and push back
      ExplicitFilteringNames nameStruct(fieldNameNode.as<std::string>(), 
                                        expFieldNameNode.as<std::string>(),
                                        fieldSizeNode.as<int>(),
                                        fieldStateSize);
      explicitFilteringNamesVec_.push_back(nameStruct);
    }
  }

  const YAML::Node y_rspecs = expect_sequence(y_explicit, "residual_specifications", true);

  // extract the sequence of types
  if (y_rspecs) {
    for (size_t ispec = 0; ispec < y_rspecs.size(); ++ispec) {
      const YAML::Node y_rspec = y_rspecs[ispec] ;
      
      // find the name, size and type
      const YAML::Node fieldNameNode = y_rspec["field_name"];
      const YAML::Node fieldSizeNode = y_rspec["field_size"];
      const YAML::Node fieldTypeNode = y_rspec["field_type"];
      
      if ( !fieldNameNode ) 
        throw std::runtime_error("ExplicitFiltering::load():error(), field_name must be provided");
      
      if ( !fieldSizeNode ) 
        throw std::runtime_error("ExplicitFiltering::load():error(), field_size must be provided");
      
      if ( fieldTypeNode ) {
        std::string fieldType = fieldTypeNode.as<std::string>();
        if ( fieldType != "node_rank" )
          throw std::runtime_error("ExplicitFiltering::load():error() only supports nodal_rank types");
      }
      
      // create the struct and push back
      ExplicitFilteringNames nameStruct(fieldNameNode.as<std::string>(), 
                                        "na",
                                        fieldSizeNode.as<int>(),
                                        1);
      residualNamesVec_.push_back(nameStruct);
    }
  }
  
  // provide output for the explicit names to be registered under STK (and, therefore, available for IO)
  if ( explicitFilteringNamesVec_.size() > 0 ) {
    NaluEnv::self().naluOutputP0() << "A number of explicitly filter quantities will be created: " << std::endl;
    NaluEnv::self().naluOutputP0() << "  " << "explicit_filter" << std::endl;
    for ( size_t k = 0; k < explicitFilteringNamesVec_.size(); ++k ) {
      NaluEnv::self().naluOutputP0() << "  " << explicitFilteringNamesVec_[k].expFieldName_ << std::endl;
    }
  }
  
  // provide output for the residual names to be registered under STK (and, therefore, available for IO)
  if ( residualNamesVec_.size() > 0 ) {
    NaluEnv::self().naluOutputP0() << "A number of residual quantities will be created: " << std::endl;
    for ( size_t k = 0; k < residualNamesVec_.size(); ++k ) {
      NaluEnv::self().naluOutputP0() << "  " << residualNamesVec_[k].fieldName_ << std::endl;
    }
  }

}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
ExplicitFiltering::setup()
{
  // objective: save off parts and register fields..
  stk::mesh::MetaData &metaData = realm_.meta_data();

  for ( size_t k = 0; k < searchTargetNames_.size(); ++k ) {
    stk::mesh::Part *thePart = metaData.get_part(searchTargetNames_[k]);
    if ( NULL != thePart )
      searchParts_.push_back(thePart);
    else
      throw std::runtime_error("ExplicitFiltering: Part is null" + searchTargetNames_[k]);
  }

  // volume for all variables
  stk::mesh::Field<double> *filter =  &(metaData.declare_field<double>(stk::topology::NODE_RANK, "explicit_filter"));
  stk::mesh::put_field_on_mesh(*filter, stk::mesh::selectUnion(searchParts_), 1, nullptr);

  // provide output for the names to be registered under STK (and, therefore, available for IO)
  for ( size_t k = 0; k < explicitFilteringNamesVec_.size(); ++k ) {
    const std::string fieldName = explicitFilteringNamesVec_[k].fieldName_;
    const std::string expFieldName = explicitFilteringNamesVec_[k].expFieldName_;
    const int fieldSize = explicitFilteringNamesVec_[k].fieldSize_;
    const int fieldStateSize = explicitFilteringNamesVec_[k].fieldStateSize_;
    stk::mesh::Field<double> *expField =  &(metaData.declare_field<double>(stk::topology::NODE_RANK, expFieldName, fieldStateSize));
    stk::mesh::put_field_on_mesh(*expField, stk::mesh::selectUnion(searchParts_), fieldSize, nullptr);
    stk::mesh::Field<double> *theField = metaData.get_field<double>(stk::topology::NODE_RANK, fieldName);
    if ( NULL != theField ) {
      // create the Np1 struct and push back
      ExplicitFilteringFields fieldStruct(theField, 
                                          expField,
                                          fieldSize);
      explicitFilteringFieldsVec_.push_back(fieldStruct);

      NaluEnv::self().naluOutputP0() << "Explicitly filter quantities with state will be created: " << expField->name() << std::endl;

      // We need to filter states as well, while later accessing 
      // the field under the general field_of_state() interface
      if ( fieldStateSize > 1 ) {
        // create the N struct and push back
        ExplicitFilteringFields fieldStructN(&(theField->field_of_state(stk::mesh::StateN)), 
                                             &(expField->field_of_state(stk::mesh::StateN)),
                                             fieldSize);
        explicitFilteringFieldsVec_.push_back(fieldStructN);

        NaluEnv::self().naluOutputP0() << "Explicitly filter quantities with state will be created: " 
                                       << (expField->field_of_state(stk::mesh::StateN)).name() << std::endl;
      }

      if ( fieldStateSize > 2 ) {        
        // create the Nm1 struct and push back
        ExplicitFilteringFields fieldStructNm1(&(theField->field_of_state(stk::mesh::StateNM1)), 
                                               &(expField->field_of_state(stk::mesh::StateNM1)),
                                               fieldSize);
        explicitFilteringFieldsVec_.push_back(fieldStructNm1);

        NaluEnv::self().naluOutputP0() << "Explicitly filter quantities with state will be created: " 
                                       << (expField->field_of_state(stk::mesh::StateNM1)).name() << std::endl;
      }      
    }
    else {
      throw std::runtime_error("ExplicitFiltering::error() no underlying field by the name of: " + fieldName);
    } 
  }

  // now residuals
  for ( size_t k = 0; k < residualNamesVec_.size(); ++k ) {
    const std::string fieldName = residualNamesVec_[k].fieldName_;
    const int fieldSize = residualNamesVec_[k].fieldSize_;
    stk::mesh::Field<double> *theField =  &(metaData.declare_field<double>(stk::topology::NODE_RANK, fieldName));
    stk::mesh::put_field_on_mesh(*theField, stk::mesh::selectUnion(searchParts_), fieldSize, nullptr);

    // create the Np1 struct and push back
    ExplicitFilteringFields fieldStruct(theField, 
                                        nullptr,
                                        fieldSize);
    residualFieldsVec_.push_back(fieldStruct);     
  }
  
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
  void
ExplicitFiltering::initialize()
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  // initialize need to ghost and elems to ghost
  needToGhostCount_ = 0;
  elemsToGhost_.clear();

  // clear the map; std::map<size_t, std::vector<stk::mesh::Entity> >
  explicitFilteringMap_.clear();

  bulkData.modification_begin();

  if ( explicitFilteringGhosting_ == NULL) {
    // create new ghosting
    std::string theGhostName = "nalu_explicit_filtering_ghosting";
    explicitFilteringGhosting_ = &bulkData.create_ghosting( theGhostName );
  }
  else {
    bulkData.destroy_ghosting(*explicitFilteringGhosting_);
  }

  bulkData.modification_end();

  // clear some of the search info
  boundingFilterVec_.clear();
  boundingElementBoxVec_.clear();
  searchKeyPair_.clear();

  // set all of the candidate elements in the search target names
  populate_candidate_elements();

  // create the ExplicitFilteringPointInfo
  create_explicit_filter_point_info_map();

  // coarse search
  determine_elems_to_ghost();

  // manage ghosting
  manage_ghosting();

  // complete filling in the set of elements connected to the centroid
  complete_search();
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ExplicitFiltering::execute()
{
  // meta/bulk data and nDim
  stk::mesh::MetaData &metaData = realm_.meta_data();
  stk::mesh::BulkData &bulkData = realm_.bulk_data();
  const int nDim = metaData.spatial_dimension();

  // extract fields
  VectorFieldType *coordinates
    = metaData.get_field<double>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // extract filter
  stk::mesh::Field<double> *explicitFilter 
    = metaData.get_field<double>(stk::topology::NODE_RANK, "explicit_filter");

  // zero assembled field(s)
  field_fill( metaData, bulkData, 0.0, *explicitFilter, realm_.get_activate_aura());  
  for ( size_t k = 0; k < explicitFilteringFieldsVec_.size(); ++k ) {
    field_fill( metaData, bulkData, 0.0, *explicitFilteringFieldsVec_[k].expField_, realm_.get_activate_aura());  
  }

  // parallel communicate data to the ghosted elements; again can communicate points to element ranks
  if ( NULL != explicitFilteringGhosting_ ) {
    std::vector< const stk::mesh::FieldBase *> ghostFieldVec;
    // fields that are needed
    ghostFieldVec.push_back(coordinates);
    // any user fields
    for ( size_t k = 0; k < explicitFilteringFieldsVec_.size(); ++k ) {
      ghostFieldVec.push_back(explicitFilteringFieldsVec_[k].theField_);
    }
    stk::mesh::communicate_field_data(*explicitFilteringGhosting_, ghostFieldVec);
  }
  
  // loop over map and assemble source terms
  std::map<stk::mesh::Entity, std::vector<stk::mesh::Entity> >::iterator iterPoint;
  for (iterPoint  = explicitFilteringMap_.begin();
       iterPoint != explicitFilteringMap_.end();
       ++iterPoint) {
    
    const stk::mesh::Entity currentNode = (*iterPoint).first;
    const std::vector<stk::mesh::Entity> &elemVec = (*iterPoint).second;

    // zero the volume
    double filterVolume = 0.0;
    for ( size_t k = 0; k < elemVec.size(); ++k ) {
      
      // extract element, nodal relationships, and nodes per element
      stk::mesh::Entity currentElem = elemVec[k];
      int nodesPerElement = bulkData.num_nodes(currentElem);

      // resize coordinates/scv and gather for subsequent volume calculation
      ws_coordinates_.resize(nodesPerElement*nDim);

      // FIXME: We may not need to gather anything else besides coordinates if we use an element-mean (no interp)
      gather_field(nDim, &ws_coordinates_[0], *coordinates, bulkData.begin_nodes(currentElem), nodesPerElement);

      // compute the volume; operates on ws_coordinates_
      double currentElemVolume = compute_volume(nDim, currentElem, bulkData);
      filterVolume += currentElemVolume;

      // loop over fields to locally assemble element average to node (might be expensive, FIXME)
      for ( size_t k = 0; k < explicitFilteringFieldsVec_.size(); ++k ) {
        const int fieldSize  = explicitFilteringFieldsVec_[k].fieldSize_;
        const stk::mesh::Field<double> *theField = explicitFilteringFieldsVec_[k].theField_;
        const stk::mesh::Field<double> *expField = explicitFilteringFieldsVec_[k].expField_;
        std::vector<double> ws_elemAveragedField(fieldSize, 0.0);
        increment_elem_mean(fieldSize, &ws_elemAveragedField[0], *theField, bulkData.begin_nodes(currentElem),
                            nodesPerElement, currentElemVolume);
        double *expF = stk::mesh::field_data(*expField, currentNode );
        for ( int j = 0; j < fieldSize; ++j )
          expF[j] += ws_elemAveragedField[j];
      }
    }
    
    // assemble explicit filter
    double *eF = stk::mesh::field_data(*explicitFilter, currentNode );
    *eF = filterVolume;
  }

  // parallel assemble; filter + user explicit filtering fields
  std::vector<const stk::mesh::FieldBase*> sumFieldVec;
  sumFieldVec.push_back(explicitFilter);
  for ( size_t k = 0; k < explicitFilteringFieldsVec_.size(); ++k ) {
    sumFieldVec.push_back(explicitFilteringFieldsVec_[k].expField_);
  }
  stk::mesh::parallel_sum(bulkData, sumFieldVec);
  
  // periodic update; filter + user fields
  if ( realm_.hasPeriodic_) {
    realm_.periodic_field_update(explicitFilter, 1);
    for ( size_t k = 0; k < explicitFilteringFieldsVec_.size(); ++k ) {
      realm_.periodic_field_update(explicitFilteringFieldsVec_[k].expField_, explicitFilteringFieldsVec_[k].fieldSize_);
    }
  }
  
  // normalize by filter 
  stk::mesh::Selector s_locally_owned = metaData.locally_owned_part()
    &stk::mesh::selectUnion(searchParts_);
  
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_locally_owned );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    
    stk::mesh::Bucket & b = **ib;
    
    const stk::mesh::Bucket::size_type length   = b.size();
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      // get node
      stk::mesh::Entity node = b[k];
      
      const double inv_eF = 1.0/(*stk::mesh::field_data(*explicitFilter, node ));
      for ( size_t k = 0; k < explicitFilteringFieldsVec_.size(); ++k ) {
        const int fieldSize  = explicitFilteringFieldsVec_[k].fieldSize_;
        const stk::mesh::Field<double> *expField = explicitFilteringFieldsVec_[k].expField_;
        double *expF = stk::mesh::field_data(*expField, node );
        for ( int j = 0; j < fieldSize; ++j ) {
          expF[j] *= inv_eF; 
        }
      }
    }
  }
  
  // parallel communicate owned to shared
  std::vector< const stk::mesh::FieldBase *> sumOwnedToSharedVec;
  for ( size_t k = 0; k < explicitFilteringFieldsVec_.size(); ++k ) {
    sumOwnedToSharedVec.push_back(explicitFilteringFieldsVec_[k].expField_);
  }
  stk::mesh::copy_owned_to_shared( bulkData, sumOwnedToSharedVec);
  
  // finally, debug output
  if ( debugOutput_ ) {
    // loop over map and provide debug
    std::map<stk::mesh::Entity, std::vector<stk::mesh::Entity> >::iterator iterPoint;
    for (iterPoint  = explicitFilteringMap_.begin();
         iterPoint != explicitFilteringMap_.end();
         ++iterPoint) {
      
      const stk::mesh::Entity currentNode = (*iterPoint).first;
      const std::vector<stk::mesh::Entity> &elemVec = (*iterPoint).second;
      
      NaluEnv::self().naluOutput() << "=====================================================" << std::endl;      
      NaluEnv::self().naluOutput() << "Current Node: " << bulkData.identifier(currentNode) 
                                   << " has element vec count: " << elemVec.size() <<  " with element ids: " << std::endl;
      for ( size_t k = 0; k < elemVec.size(); ++k ) {
        NaluEnv::self().naluOutput() << bulkData.identifier(elemVec[k]) << " "; 
      }
      NaluEnv::self().naluOutput() << std::endl;
    }
  }
  
  // residual-based calculation is currently in a very specific state
  if ( residualNamesVec_.size() > 0 ) {
    // compute the non-filtered residuals
    VectorFieldType *velocity = metaData.get_field<double>(stk::topology::NODE_RANK, "velocity");
    ScalarFieldType *density = metaData.get_field<double>(stk::topology::NODE_RANK, "density");
    ScalarFieldType *viscosity = metaData.get_field<double>(stk::topology::NODE_RANK, "viscosity");
    VectorFieldType *residual = metaData.get_field<double>(stk::topology::NODE_RANK, "residual");
    
    ScalarFieldType *pressure = metaData.get_field<double>(stk::topology::NODE_RANK, "pressure");
    VectorFieldType *dudx = metaData.get_field<double>(stk::topology::NODE_RANK, "dudx");
    compute_scv_residual(&(velocity->field_of_state(stk::mesh::StateNP1)),
                         &(velocity->field_of_state(stk::mesh::StateN)),
                         &(velocity->field_of_state(stk::mesh::StateNM1)),
                         &(density->field_of_state(stk::mesh::StateNP1)),
                         &(density->field_of_state(stk::mesh::StateN)),
                         &(density->field_of_state(stk::mesh::StateNM1)),
                         viscosity,
                         pressure,
                         dudx,
                         coordinates,
                         residual);
  }
  
  if ( residualNamesVec_.size() > 1 ) {
    // compute the filtered residuals
    VectorFieldType *evelocity = metaData.get_field<double>(stk::topology::NODE_RANK, "explicit_velocity");
    ScalarFieldType *edensity = metaData.get_field<double>(stk::topology::NODE_RANK, "explicit_density");
    ScalarFieldType *eviscosity = metaData.get_field<double>(stk::topology::NODE_RANK, "explicit_viscosity");
    VectorFieldType *eresidual = metaData.get_field<double>(stk::topology::NODE_RANK, "explicit_residual");    
    ScalarFieldType *epressure = metaData.get_field<double>(stk::topology::NODE_RANK, "explicit_pressure");
    VectorFieldType *edudx = metaData.get_field<double>(stk::topology::NODE_RANK, "explicit_dudx");
    compute_scv_residual(&(evelocity->field_of_state(stk::mesh::StateNP1)),
                         &(evelocity->field_of_state(stk::mesh::StateN)),
                         &(evelocity->field_of_state(stk::mesh::StateNM1)),
                         &(edensity->field_of_state(stk::mesh::StateNP1)),
                         &(edensity->field_of_state(stk::mesh::StateN)),
                         &(edensity->field_of_state(stk::mesh::StateNM1)),
                         eviscosity,
                         epressure,
                         edudx,
                         coordinates,
                         eresidual);
    
  }
  
  if ( residualNamesVec_.size() > 2 ) {
    throw std::runtime_error("ExplicitFiltering::execute() only two entries supported");
  }

}
  
//--------------------------------------------------------------------------
//-------- populate_candidate_elements -------------------------------------
//--------------------------------------------------------------------------
void
ExplicitFiltering::populate_candidate_elements()
{
  stk::mesh::MetaData &metaData = realm_.meta_data();
  stk::mesh::BulkData &bulkData = realm_.bulk_data();

  const int nDim = metaData.spatial_dimension();

  // fields
  VectorFieldType *coordinates = metaData.get_field<double>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // point data structures
  Point minCorner, maxCorner;

  // selector and bucket loop
  stk::mesh::Selector s_locally_owned = metaData.locally_owned_part()
    &stk::mesh::selectUnion(searchParts_);

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned );

  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {

    stk::mesh::Bucket & b = **ib;

    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get element
      stk::mesh::Entity elem = b[k];

      // initialize max and min
      for (int j = 0; j < nDim; ++j ) {
        minCorner[j] = +1.0e16;
        maxCorner[j] = -1.0e16;
      }

      // extract elem_node_relations
      stk::mesh::Entity const* elem_node_rels = bulkData.begin_nodes(elem);
      const int num_nodes = bulkData.num_nodes(elem);

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];

        // pointers to real data
        const double * coords = stk::mesh::field_data(*coordinates, node );

        // check max/min
        for ( int j = 0; j < nDim; ++j ) {
          minCorner[j] = std::min(minCorner[j], coords[j]);
          maxCorner[j] = std::max(maxCorner[j], coords[j]);
        }
      }

      // setup ident
      stk::search::IdentProc<uint64_t,int> theIdent(bulkData.identifier(elem), NaluEnv::self().parallel_rank());

      // create the bounding point box and push back
      boundingElementBox theBox(Box(minCorner,maxCorner), theIdent);
      boundingElementBoxVec_.push_back(theBox);
    }
  }
}

//--------------------------------------------------------------------------
//-------- determine_elems_to_ghost ----------------------------------------
//--------------------------------------------------------------------------
void
ExplicitFiltering::determine_elems_to_ghost()
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  stk::search::coarse_search(boundingFilterVec_, boundingElementBoxVec_, searchMethod_,
                             NaluEnv::self().parallel_comm(), searchKeyPair_);

  // sort to avoid possible elemsToGhost_ ordering in change_ghosting() or actuatorLinePointInfo->nodeVec_.insert(node)
  std::sort(searchKeyPair_.begin(), searchKeyPair_.end());
  
  // lowest effort is to ghost elements to the owning rank of the point; can just as easily do the opposite
  std::vector<std::pair<boundingElementBox::second_type, boundingElementBox::second_type> >::const_iterator ii;
  for( ii=searchKeyPair_.begin(); ii!=searchKeyPair_.end(); ++ii ) {

    const uint64_t theBox = ii->second.id();
    unsigned theRank = NaluEnv::self().parallel_rank();
    const unsigned pt_proc = ii->first.proc();
    const unsigned box_proc = ii->second.proc();
    
    if ( (box_proc == theRank) && (pt_proc != theRank) ) {

      // Send box to pt proc

      // find the element
      stk::mesh::Entity theElemMeshObj = bulkData.get_entity(stk::topology::ELEMENT_RANK, theBox);
      if ( !(bulkData.is_valid(theElemMeshObj)) )
        throw std::runtime_error("no valid entry for element");

      // new element to ghost counter
      needToGhostCount_++;

      // deal with elements to push back to be ghosted
      stk::mesh::EntityProc theElemPair(theElemMeshObj, pt_proc);
      elemsToGhost_.push_back(theElemPair);
    }
  }
}

//--------------------------------------------------------------------------
//-------- create_explicit_filter_point_info_map -----------------------------
//--------------------------------------------------------------------------
void
ExplicitFiltering::create_explicit_filter_point_info_map() {

  stk::mesh::MetaData &metaData = realm_.meta_data();
  stk::mesh::BulkData &bulkData = realm_.bulk_data();

  const int nDim = metaData.spatial_dimension();

  // abstract out filterSize to avoid logic in a deep loop
  double filterSizeAbstract[3] = {0.0,0.0,0.0};
  filterSizeAbstract[0] = filterSize_.x_;
  filterSizeAbstract[1] = filterSize_.y_;
  if ( nDim > 2 )
    filterSizeAbstract[2] = filterSize_.z_;
  
  // fields
  VectorFieldType *coordinates = metaData.get_field<double>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // define a point that will hold the min/max
  Point minCorner, maxCorner;

  // selector and bucket loop
  stk::mesh::Selector s_locally_owned = metaData.locally_owned_part()
    &stk::mesh::selectUnion(searchParts_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_locally_owned );

  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {

    stk::mesh::Bucket & b = **ib;

    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get node
      stk::mesh::Entity node = b[k];

      // pointers to real data
      const double * coords = stk::mesh::field_data(*coordinates, node );

      // determine min/max locations; center a box about the nodal coordinates
      for ( int j = 0; j < nDim; ++j ) {
        const double dxj = 0.5*filterSizeAbstract[j];
        minCorner[j] = coords[j] - dxj;
        maxCorner[j] = coords[j] + dxj;
      }
      
      if ( debugOutput_ ) {
        NaluEnv::self().naluOutput() << "Current Node: " << bulkData.identifier(node) << " pt coords and min/max corner:" << std::endl;
        for ( int j = 0; j < nDim; ++j ) {
          NaluEnv::self().naluOutput() << coords[j] << std::endl;
          NaluEnv::self().naluOutput() << minCorner[j] << "/" << maxCorner[j] << std::endl;
        }
      }

      // setup ident
      stk::search::IdentProc<uint64_t,int> theIdent(bulkData.identifier(node), NaluEnv::self().parallel_rank());
      
      // create the bounding box for the filter and push back
      boundingElementBox theBox(Box(minCorner,maxCorner), theIdent);
      boundingFilterVec_.push_back(theBox);

      // create the map with an empty vector (to be filled in complete_search)
      std::vector<stk::mesh::Entity> elemVec;
      explicitFilteringMap_[node] = elemVec;
    }
  }
}
  
//--------------------------------------------------------------------------
//-------- manage_ghosting -------------------------------------------------
//--------------------------------------------------------------------------
void
ExplicitFiltering::manage_ghosting()
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  // check for ghosting need
  uint64_t g_needToGhostCount = 0;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &needToGhostCount_, &g_needToGhostCount, 1);
  if (g_needToGhostCount > 0) {
    NaluEnv::self().naluOutputP0() << "ExplicitFiltering alg will ghost a number of entities: "
                                   << g_needToGhostCount  << std::endl;    
    bulkData.modification_begin();
    bulkData.change_ghosting( *explicitFilteringGhosting_, elemsToGhost_);
    bulkData.modification_end();
  }
  else {
    NaluEnv::self().naluOutputP0() << "ExplicitFiltering alg will NOT ghost entities: " << std::endl;
  }
}

//--------------------------------------------------------------------------
//-------- complete_search -------------------------------------------------
//--------------------------------------------------------------------------
void
ExplicitFiltering::complete_search()
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  // now proceed with the standard search
  std::vector<std::pair<boundingElementBox::second_type, boundingElementBox::second_type> >::const_iterator ii;
  for( ii=searchKeyPair_.begin(); ii!=searchKeyPair_.end(); ++ii ) {

    const uint64_t thePt = ii->first.id();
    const uint64_t theBox = ii->second.id();
    const unsigned theRank = NaluEnv::self().parallel_rank();
    const unsigned pt_proc = ii->first.proc();

    // check if I own the point...
    if ( theRank == pt_proc ) {

      // yes, I own the point...

      // proceed as required; all elements should have already been ghosted via the coarse search
      stk::mesh::Entity elem = bulkData.get_entity(stk::topology::ELEMENT_RANK, theBox);
      if ( !(bulkData.is_valid(elem)) )
        throw std::runtime_error("no valid entry for element");

      // find the point data structure
      stk::mesh::Entity node = bulkData.get_entity(stk::topology::NODE_RANK, thePt);
      if ( !(bulkData.is_valid(node)) )
        throw std::runtime_error("no valid entry for node");

      // extract the point object and push back the element
      std::vector<stk::mesh::Entity> &elemVec = explicitFilteringMap_[node];
      elemVec.push_back(elem);
    }
    else {
      // not this proc's issue
    }
  }

  if ( debugOutput_ ) {
    // global sum on local count
    size_t l_count = explicitFilteringMap_.size();
    size_t g_count = 0;
    stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
    stk::all_reduce_sum(comm, &l_count, &g_count, 1);

    // output
    NaluEnv::self().naluOutputP0() << "Total size of Map: " << g_count << std::endl;
    
    // now provide total node count
    l_count = 0;
    g_count = 0;

    stk::mesh::MetaData &metaData = realm_.meta_data();

    stk::mesh::Selector s_locally_owned = metaData.locally_owned_part()
      &stk::mesh::selectUnion(searchParts_);
    
    stk::mesh::BucketVector const& node_buckets =
      realm_.get_buckets( stk::topology::NODE_RANK, s_locally_owned );
    
    for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
          ib != node_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib;      
      l_count += b.size();
    }

    // output
    stk::all_reduce_sum(comm, &l_count, &g_count, 1);
    NaluEnv::self().naluOutputP0() << "Global node count: " << g_count << std::endl;
  }
  
}

//--------------------------------------------------------------------------
//-------- gather_field ----------------------------------------------------
//--------------------------------------------------------------------------
void
ExplicitFiltering::gather_field(
  const int &sizeOfField,
  double *fieldToFill,
  const stk::mesh::FieldBase &stkField,
  stk::mesh::Entity const* elem_node_rels,
  const int &nodesPerElement)
{
  for ( int ni = 0; ni < nodesPerElement; ++ni ) {
    stk::mesh::Entity node = elem_node_rels[ni];
    const double * theField = (double*)stk::mesh::field_data(stkField, node );
    for ( int j = 0; j < sizeOfField; ++j ) {
      const int offSet = ni*sizeOfField+j;
      fieldToFill[offSet] = theField[j];
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_volume --------------------------------------------------
//--------------------------------------------------------------------------
double
ExplicitFiltering::compute_volume(
  const int &nDim,
  stk::mesh::Entity elem,
  const stk::mesh::BulkData & bulkData)
{
  // extract master element from the bucket in which the element resides
  const stk::topology &elemTopo = bulkData.bucket(elem).topology();
  MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(elemTopo);
  const int numScvIp = meSCV->numIntPoints_;

  // compute scv for this element
  ws_scv_volume_.resize(numScvIp);
  double scv_error = 0.0;
  meSCV->determinant(1, &ws_coordinates_[0], &ws_scv_volume_[0], &scv_error);
  
  double elemVolume = 0.0;
  for ( int ip = 0; ip < numScvIp; ++ip ) {
    elemVolume += ws_scv_volume_[ip];
  }
  return elemVolume;
}

//--------------------------------------------------------------------------
//-------- increment_elem_mean ---------------------------------------------
//--------------------------------------------------------------------------
void
ExplicitFiltering::increment_elem_mean(
  const int &sizeOfField,
  double *elemAveragedField,
  const stk::mesh::FieldBase &stkField,
  stk::mesh::Entity const* elem_node_rels,
  const int &nodesPerElement,
  const double &elemVolume)
{
  // FIXME: Should probably gather the nodal field, interpolate to scv and scale by scv_vol
  // for now, we are computing the simple mean (sum_field/npe) over the element and accumulating
  const double invScaling = 1.0/(double)nodesPerElement;
  for ( int ni = 0; ni < nodesPerElement; ++ni ) {
    stk::mesh::Entity node = elem_node_rels[ni];
    const double * theField = (double*)stk::mesh::field_data(stkField, node );
    for ( int j = 0; j < sizeOfField; ++j ) {
      elemAveragedField[j] += theField[j]*invScaling*elemVolume;
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_scv_residual --------------------------------------------
//--------------------------------------------------------------------------
void
ExplicitFiltering::compute_scv_residual(
  const stk::mesh::FieldBase *velocityNp1_,
  const stk::mesh::FieldBase *velocityN_,
  const stk::mesh::FieldBase *velocityNm1_,
  const stk::mesh::FieldBase *densityNp1_,
  const stk::mesh::FieldBase *densityN_,
  const stk::mesh::FieldBase *densityNm1_,
  const stk::mesh::FieldBase *viscosity_,
  const stk::mesh::FieldBase *pressure_,
  const stk::mesh::FieldBase *Gjui_,
  const stk::mesh::FieldBase *coordinates_,
  stk::mesh::FieldBase *residual_)
{
  NaluEnv::self().naluOutputP0() << "Nalu residual acting on: " << residual_->name() << " using: " << std::endl;
  NaluEnv::self().naluOutputP0() << velocityNp1_->name() << " " 
                                 << velocityN_->name() << " "
                                 << velocityNm1_->name() << " "
                                 << densityNp1_->name() << " "
                                 << densityN_->name() << " "
                                 << densityNm1_->name() << " "
                                 << viscosity_->name() << " "
                                 << pressure_->name() << " " 
                                 << Gjui_->name() << std::endl;
  
  // compute assembled residual to nodes for momentum system 
  stk::mesh::MetaData &metaData = realm_.meta_data();
  stk::mesh::BulkData &bulkData = realm_.bulk_data();
  const int nDim = metaData.spatial_dimension();

  // extract constants
  const double rhoRef = realm_.solutionOptions_->referenceDensity_;
  const std::array<double,3>& gravity = realm_.solutionOptions_->get_gravity_vector();
  const double gamma1 = realm_.get_gamma1();
  const double gamma2 = realm_.get_gamma2();
  const double gamma3 = realm_.get_gamma3();
  const double dt = realm_.get_time_step();
  const double includeDivU = realm_.get_divU();
  const double twoThirds = 2.0/3.0;
  
  // subtract out continuity factor and kd
  const double ncFac = 1.0;
  double w_kd[3][3] = {{1.0, 0.0, 0.0}, 
                       {0.0, 1.0, 0.0}, 
                       {0.0, 0.0, 1.0}};

  // zero residual
  field_fill( metaData, bulkData, 0.0, *residual_, realm_.get_activate_aura());  

  std::vector<double> ws_residual;
  std::vector<double> ws_velocityNp1;
  std::vector<double> ws_velocityN;
  std::vector<double> ws_velocityNm1;
  std::vector<double> ws_coordinates;
  std::vector<double> ws_densityNp1;
  std::vector<double> ws_densityN;
  std::vector<double> ws_densityNm1;
  std::vector<double> ws_viscosity;
  std::vector<double> ws_pressure;
  std::vector<double> ws_Gjui;

  // geometry related to populate
  std::vector<double> ws_dndx;
  std::vector<double> ws_deriv;
  std::vector<double> ws_det_j;
  std::vector<double> ws_shape_function;
  std::vector<double> ws_scvol;

  // ip values
  std::vector<double> uNp1Ip(nDim);
  std::vector<double> uNIp(nDim);
  std::vector<double> uNm1Ip(nDim);
  std::vector<double> dpdxIp(nDim);
  std::vector<double> dFdxAdv(nDim);
  std::vector<double> dFdxDiff(nDim);

  // pointers for fast access
  double *p_uNp1Ip = &uNp1Ip[0];
  double *p_uNIp = &uNIp[0];
  double *p_uNm1Ip = &uNm1Ip[0];
  double *p_dpdxIp = &dpdxIp[0];
  double *p_dFdxAdv = &dFdxAdv[0];
  double *p_dFdxDiff = &dFdxDiff[0];
  
  stk::mesh::Selector s_locally_owned_union = metaData.locally_owned_part()
    & stk::mesh::selectUnion(searchParts_) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meSCV->nodesPerElement_;
    const int numScvIp = meSCV->numIntPoints_;
    const int *ipNodeMap = meSCV->ipNodeMap();

    // algorithm related
    ws_residual.resize(nodesPerElement*nDim);
    ws_velocityNp1.resize(nodesPerElement*nDim);
    ws_velocityN.resize(nodesPerElement*nDim);
    ws_velocityNm1.resize(nodesPerElement*nDim);
    ws_coordinates.resize(nodesPerElement*nDim);
    ws_Gjui.resize(nDim*nDim*nodesPerElement);
    ws_pressure.resize(nodesPerElement);
    ws_densityNp1.resize(nodesPerElement);
    ws_densityN.resize(nodesPerElement);
    ws_densityNm1.resize(nodesPerElement);
    ws_viscosity.resize(nodesPerElement);
    ws_dndx.resize(nDim*numScvIp*nodesPerElement);
    ws_deriv.resize(nDim*numScvIp*nodesPerElement);
    ws_det_j.resize(numScvIp);
    ws_shape_function.resize(numScvIp*nodesPerElement);

    ws_scvol.resize(numScvIp);

    // pointer to residual
    double *p_residual = &ws_residual[0];

    double *p_velocityNp1 = &ws_velocityNp1[0];
    double *p_velocityN = &ws_velocityN[0];
    double *p_velocityNm1 = &ws_velocityNm1[0];
    double *p_coordinates = &ws_coordinates[0];
    double *p_Gjui = &ws_Gjui[0];
    double *p_pressure = &ws_pressure[0];
    double *p_densityNp1 = &ws_densityNp1[0];
    double *p_densityN = &ws_densityN[0];
    double *p_densityNm1 = &ws_densityNm1[0];
    double *p_viscosity = &ws_viscosity[0];
    double *p_dndx = &ws_dndx[0];
    double *p_shape_function = &ws_shape_function[0];

    double *p_scvol = &ws_scvol[0];
    
    // extract shape function
    meSCV->shape_fcn(&p_shape_function[0]);
       
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const * node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);

      // sanity check on num nodes
      STK_ThrowAssert( num_nodes == nodesPerElement );

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // pointers to real data
        const double * uNp1     =  (double*)stk::mesh::field_data(*velocityNp1_, node);
        const double * uN       =  (double*)stk::mesh::field_data(*velocityN_, node);
        const double * uNm1     =  (double*)stk::mesh::field_data(*velocityNm1_, node);
        const double * coords   =  (double*)stk::mesh::field_data(*coordinates_, node);
        const double * Gjui     =  (double*)stk::mesh::field_data(*Gjui_, node);
        const double pressure   = *(double*)stk::mesh::field_data(*pressure_, node);
        const double rhoNp1     = *((double*)stk::mesh::field_data(*densityNp1_, node));
        const double rhoN       = *((double*)stk::mesh::field_data(*densityN_, node));
        const double rhoNm1     = *((double*)stk::mesh::field_data(*densityNm1_, node));
        const double mu         = *((double*)stk::mesh::field_data(*viscosity_, node));

        // gather scalars
        p_densityNp1[ni] = rhoNp1;
        p_densityN[ni] = rhoN;
        p_densityNm1[ni] = rhoNm1;
        p_viscosity[ni] = mu;
        p_pressure[ni] = pressure;

        // gather vectors
        const int niNdim = ni*nDim;
        const int niNdimNdim = niNdim*nDim;
        for ( int i = 0; i < nDim; ++i ) {
          p_velocityNp1[niNdim+i] = uNp1[i];
          p_velocityN[niNdim+i] = uN[i];
          p_velocityNm1[niNdim+i] = uNm1[i];
          p_coordinates[niNdim+i] = coords[i];
          // gather tensor
          const int iNdim = i*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            p_Gjui[niNdimNdim+iNdim+j] = Gjui[iNdim+j];
          }
        }
      }
      
      // compute geometry
      double scv_error = 0.0;
      meSCV->determinant(1, &p_coordinates[0], &p_scvol[0], &scv_error);

      // compute dndx
      meSCV->grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_det_j[0], &scv_error);

      // zero residual for this element; residual follows matrix assembly rules
      for ( int j = 0; j < nodesPerElement*nDim; ++j )
        p_residual[j] = 0.0;
      
      // loop over SCV only; time, advection, duffusion, pressure gradient, and buoyancy accumulation
      for ( int ip = 0; ip < numScvIp; ++ip ) {

        // save off some offsets
        const int nn = ipNodeMap[ip];
        const int ipNpe = ip*nodesPerElement;

        // zero everything that we will need at the IP
        double rhoNp1Ip = 0.0;
        double rhoNIp = 0.0;
        double rhoNm1Ip = 0.0;
        double divU = 0.0;
        double dFdxCont = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          p_uNp1Ip[j] = 0.0;
          p_uNIp[j] = 0.0;
          p_uNm1Ip[j] = 0.0;
          p_dFdxAdv[j] = 0.0;
          p_dFdxDiff[j] = 0.0;
          p_dpdxIp[j] = 0.0;
        }

        // compute everything other than diffusion
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {  
          const double r = p_shape_function[ipNpe+ic];        
          
          rhoNp1Ip += r*p_densityNp1[ic];
          rhoNIp += r*p_densityN[ic];
          rhoNm1Ip += r*p_densityNm1[ic];

          const double pIc = p_pressure[ic];
          const double rhoNp1Ic = p_densityNp1[ic];
          const int icNdim = ic*nDim;
          const int offSetDnDx = nDim*nodesPerElement*ip + icNdim;
          for ( int j = 0; j < nDim; ++j ) {
            // interpolation
            const double uNp1j =p_velocityNp1[icNdim+j]; 
            p_uNp1Ip[j] += r*uNp1j;
            p_uNIp[j] += r*p_velocityN[icNdim+j];
            p_uNm1Ip[j] += r*p_velocityNm1[icNdim+j];
            // derivatives
            const double dnj = p_dndx[offSetDnDx+j];
            p_dpdxIp[j] += pIc*dnj;
            divU += uNp1j*dnj;
            dFdxCont += rhoNp1Ic*uNp1j*dnj;
          }
        }
        
        // full continuity residual
        const double contRes = (gamma1*rhoNp1Ip + gamma2*rhoNIp + gamma3*rhoNm1Ip)/dt + dFdxCont;
        
        // advection and diffusion terms
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const int icNdim = ic*nDim;
          const int icNdimNdim = icNdim*nDim;
          const int offSetDnDx = nDim*nodesPerElement*ip + icNdim;
          const double rhoNp1Ic = p_densityNp1[ic];
          const double muIc = p_viscosity[ic];
          for ( int i = 0; i < nDim; ++i ) {
            for ( int j = 0; j < nDim; ++j ) {
              const double dnj = p_dndx[offSetDnDx+j];
              p_dFdxAdv[i] += rhoNp1Ic*p_velocityNp1[icNdim+j]*p_velocityNp1[icNdim+i]*dnj;
              p_dFdxDiff[i] += muIc*(p_Gjui[icNdimNdim+i*nDim+j] + p_Gjui[icNdimNdim+j*nDim+i] 
                                     - twoThirds*divU*w_kd[i][j]*includeDivU)*dnj;
            }
          }
        }
        
        // offset
        const int nnNdim = nn*nDim;
        
        // save off scvol
        const double scvIp = p_scvol[nn];

        for ( int i = 0; i < nDim; ++i ) {
          const int indexNN = nnNdim + i;
          const double time = (gamma1*rhoNp1Ip*p_uNp1Ip[i] 
                               + gamma2*rhoNIp*p_uNIp[i] 
                               + gamma3*rhoNm1Ip*p_uNm1Ip[i])/dt;
          const double buoyancy = (rhoNp1Ip-rhoRef)*gravity[i];
          p_residual[indexNN] += 
            (time + p_dFdxAdv[i] - p_dFdxDiff[i] + p_dpdxIp[i] - buoyancy - ncFac*contRes*p_uNp1Ip[i])*scvIp;
        }
      }
       
      // scatter residual
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // off set
        const int niNdim = ni*nDim;
        
        // pointer to the node
        double * residual =  (double*)stk::mesh::field_data(*residual_, node);
        for ( int i = 0; i < nDim; ++i ) {
          residual[i] -= p_residual[niNdim+i];
        }
      } 
    }
  }
  
  // parallel assemble; first
  std::vector<const stk::mesh::FieldBase*> sumFieldVecFirst;
  sumFieldVecFirst.push_back(residual_);
  stk::mesh::parallel_sum(bulkData, sumFieldVecFirst);
  
  // periodic sum
  if ( realm_.hasPeriodic_) {
    realm_.periodic_field_update(residual_, nDim);
  }

  // choice to normalize
  if ( normalizeResidual_ ) {
    ScalarFieldType *dualNodalVolume = metaData.get_field<double>(stk::topology::NODE_RANK, "dual_nodal_volume");
    field_normalize(metaData, bulkData, *dualNodalVolume, *residual_, realm_.get_activate_aura());
  }
}

} // namespace nalu
} // namespace Sierra
