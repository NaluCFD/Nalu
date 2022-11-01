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
    debugOutput_(false)
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
  }

  // extract the sequence of types
  const YAML::Node y_specs = expect_sequence(y_explicit, "specifications", true);
  if (y_specs) {
    for (size_t ispec = 0; ispec < y_specs.size(); ++ispec) {
      const YAML::Node y_spec = y_specs[ispec] ;
      
      // find the name, size and type
      const YAML::Node fieldNameNode = y_spec["field_name"];
      const YAML::Node expFieldNameNode = y_spec["explicit_field_name"];
      const YAML::Node fieldSizeNode = y_spec["field_size"];
      const YAML::Node fieldTypeNode = y_spec["field_type"];
      
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

      // create the struct and push back
      ExplicitFilteringNames nameStruct(fieldNameNode.as<std::string>(), 
                                        expFieldNameNode.as<std::string>(),
                                        fieldSizeNode.as<int>());
      explicitFilteringNamesVec_.push_back(nameStruct);
    }
  }

  // provide output for the names to be registered under STK (and, therefore, available for IO)
  NaluEnv::self().naluOutputP0() << "A number of explicitly filter quantities will be created: " << std::endl;
  NaluEnv::self().naluOutputP0() << "  " << "explicit_filter" << std::endl;
  for ( size_t k = 0; k < explicitFilteringNamesVec_.size(); ++k ) {
    NaluEnv::self().naluOutputP0() << "  " << explicitFilteringNamesVec_[k].expFieldName_ << std::endl;
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
  stk::mesh::Field<double> *filter =  &(metaData.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "explicit_filter"));
  stk::mesh::put_field_on_mesh(*filter, stk::mesh::selectUnion(searchParts_), 1, nullptr);

  // provide output for the names to be registered under STK (and, therefore, available for IO)
  for ( size_t k = 0; k < explicitFilteringNamesVec_.size(); ++k ) {
    const std::string fieldName = explicitFilteringNamesVec_[k].fieldName_;
    const std::string expFieldName = explicitFilteringNamesVec_[k].expFieldName_;
    const int fieldSize = explicitFilteringNamesVec_[k].fieldSize_;
    stk::mesh::Field<double> *expField =  &(metaData.declare_field<double>(stk::topology::NODE_RANK, expFieldName));
    stk::mesh::put_field_on_mesh(*expField, stk::mesh::selectUnion(searchParts_), fieldSize, nullptr);
    stk::mesh::Field<double> *theField = metaData.get_field<double>(stk::topology::NODE_RANK, fieldName);
    if ( NULL != theField ) {
      // create the struct and push back
      ExplicitFilteringFields fieldStruct(theField, 
                                          expField,
                                          fieldSize);
      explicitFilteringFieldsVec_.push_back(fieldStruct);
    }
    else {
      throw std::runtime_error("ExplicitFiltering::error() no underlying field by the name of: " + fieldName);
    } 
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
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const int nDim = metaData.spatial_dimension();

  // extract fields
  VectorFieldType *coordinates
    = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // extract filter
  stk::mesh::Field<double> *explicitFilter 
    = metaData.get_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "explicit_filter");

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
  
}

//--------------------------------------------------------------------------
//-------- populate_candidate_elements -------------------------------------
//--------------------------------------------------------------------------
void
ExplicitFiltering::populate_candidate_elements()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  const int nDim = metaData.spatial_dimension();

  // fields
  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

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

  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  const int nDim = metaData.spatial_dimension();

  // abstract out filterSize to avoid logic in a deep loop
  double filterSizeAbstract[3] = {0.0,0.0,0.0};
  filterSizeAbstract[0] = filterSize_.x_;
  filterSizeAbstract[1] = filterSize_.y_;
  if ( nDim > 2 )
    filterSizeAbstract[2] = filterSize_.z_;
  
  // fields
  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

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


} // namespace nalu
} // namespace Sierra
