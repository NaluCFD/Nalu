/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ContactInfo.h>
#include <ContactManager.h>
#include <HaloInfo.h>
#include <master_element/MasterElement.h>
#include <Realm.h>
#include <NaluEnv.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/CPUTime.hpp>

// stk_search
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>

// vector and pair
#include <vector>

namespace sierra{
namespace nalu{

class HaloInfo;

//==========================================================================
// Class Definition
//==========================================================================
// Acon_ContactInfo - contains virtual edge data
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ContactInfo::ContactInfo(
   Realm &realm,
   const std::string &name,
   const double maxSearchRadius,
   const double minSearchRadius,
   const std::vector<std::string> &contactSearchBlockName,
   const double expandBoxPercentage,
   const stk::mesh::Part *contactSurfacePart,
   const std::string &searchMethodName,
   const bool clipIsoParametricCoords,
   const bool useHermiteInterpolation)
  : realm_(realm ),
    name_(name ),
    maxSearchRadius_(maxSearchRadius),
    minSearchRadius_(minSearchRadius),
    contactSurfacePart_(contactSurfacePart),
    meSCS_(NULL),
    expandBoxPercentage_(expandBoxPercentage),
    meshMotion_(realm_.has_mesh_motion()),
    searchMethod_(stk::search::BOOST_RTREE),
    clipIsoParametricCoords_(clipIsoParametricCoords),
    useHermiteInterpolation_(useHermiteInterpolation)
{
  // determine search method for this pair
  if ( searchMethodName == "boost_rtree" )
    searchMethod_ = stk::search::BOOST_RTREE;
  else if ( searchMethodName == "stk_octree" )
    searchMethod_ = stk::search::OCTREE;
  else
    NaluEnv::self().naluOutputP0() << "ContactInfo::search method not declared; will use BOOST_RTREE" << std::endl;

  // save off vector of parts for the search blocks
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  for ( size_t k = 0; k < contactSearchBlockName.size(); ++k ) {
     stk::mesh::Part *blockSearchPart = meta_data.get_part(contactSearchBlockName[k]);
     if ( NULL == blockSearchPart ) {
       throw std::runtime_error("Sorry, issue with search block specification");
     }
     contactSearchBlock_.push_back(blockSearchPart);
   }

}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ContactInfo::~ContactInfo()
{
  // delete halo info objects:
  std::map<uint64_t, HaloInfo*>::iterator ii;
  for( ii=haloInfoMap_.begin(); ii!=haloInfoMap_.end(); ++ii )
    delete (*ii).second;
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
ContactInfo::initialize()
{

  // clear some of the search info
  boundingPointVec_.clear();
  boundingElementBoxVec_.clear();
  searchKeyPair_.clear();

  // delete haloInfoMap_
  std::map<uint64_t, HaloInfo *>::iterator ii;
  for( ii=haloInfoMap_.begin(); ii!=haloInfoMap_.end(); ++ii )
    delete (*ii).second;
  haloInfoMap_.clear();

  construct_halo_state();

  if ( meshMotion_ )
    populate_halo_mesh_velocity();

  find_possible_elements();

  set_best_x();

  determine_elems_to_ghost();

}

//--------------------------------------------------------------------------
//-------- construct_halo_state --------------------------------------------
//--------------------------------------------------------------------------
void
ContactInfo::construct_halo_state()
{

  // hack due to 3d
  Point localHaloCoords;

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  // fields
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  VectorFieldType *haloAxj = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "halo_axj");
  VectorFieldType *haloDxj = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "halo_dxj");

  const int nDim = meta_data.spatial_dimension();

  stk::mesh::Selector s_locally_owned = meta_data.locally_owned_part()
    &stk::mesh::Selector(*contactSurfacePart_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_locally_owned );

  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {

    stk::mesh::Bucket & b = **ib;

    const stk::mesh::Bucket::size_type length   = b.size();

    // point to data
    const double * coords = stk::mesh::field_data(*coordinates, b);
    double * hAxj = stk::mesh::field_data(*haloAxj, b);
    double * hDxj = stk::mesh::field_data(*haloDxj, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get face
      stk::mesh::Entity node = b[k];

      HaloInfo *haloInfo = new HaloInfo(node, nDim);

      // setup ident; do something about processor count...
      stk::search::IdentProc<uint64_t,int> theIdent(bulk_data.identifier(node), NaluEnv::self().parallel_rank());
      haloInfoMap_[bulk_data.identifier(node)] = haloInfo;

      // define offset for all nodal fields that are of nDim
      const size_t offSet = k*nDim;

      // populate haloNodalCoords
      double dS = 0.0;
      for (int j = 0; j < nDim; ++j ) {
        const double xj = coords[offSet+j];
        const double dxj = hDxj[offSet+j];
        dS += dxj*dxj;
        const double hncj = xj + dxj;
        // store halo nodal coords to local to provide to bounding point
        localHaloCoords[j] = hncj;
        // save off all
        const double axj = hAxj[offSet+j];
        haloInfo->haloEdgeAreaVec_[j] = axj;
        haloInfo->nodalCoords_[j] = xj;
        haloInfo->haloNodalCoords_[j] = hncj;
      }

      haloInfo->haloEdgeDs_ = std::sqrt(dS);
      // create the bounding point box and push back
      boundingPoint thePt(localHaloCoords, theIdent);
      boundingPointVec_.push_back(thePt);
    }
  }
}

//--------------------------------------------------------------------------
//-------- populate_halo_mesh_velocity -------------------------------------
//--------------------------------------------------------------------------
void
ContactInfo::populate_halo_mesh_velocity()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // fields
  ScalarFieldType *omegaField = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "omega");

  // iterate halo face nodes
  std::map<uint64_t, HaloInfo *>::iterator iterHalo;
  for (iterHalo  = haloInfoMap_.begin();
       iterHalo != haloInfoMap_.end();
       ++iterHalo) {

    // halo info object of interest
    HaloInfo * infoObject = (*iterHalo).second;

    // extract coordinates for halo node
    const double cX = infoObject->haloNodalCoords_[0];
    const double cY = infoObject->haloNodalCoords_[1];

    // extract omega
    const double omega = *stk::mesh::field_data(*omegaField, infoObject->faceNode_);

    // set mesh velocity
    infoObject->haloMeshVelocity_[0] = -omega*cY;
    infoObject->haloMeshVelocity_[1] = +omega*cX;
  }
}

//--------------------------------------------------------------------------
//-------- determine_elems_to_ghost ----------------------------------------
//--------------------------------------------------------------------------
void
ContactInfo::determine_elems_to_ghost()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  stk::search::coarse_search(boundingPointVec_, boundingElementBoxVec_, searchMethod_, NaluEnv::self().parallel_comm(), searchKeyPair_);

  std::vector<std::pair<boundingPoint::second_type, boundingElementBox::second_type> >::const_iterator ii;
  for( ii=searchKeyPair_.begin(); ii!=searchKeyPair_.end(); ++ii ) {

    const uint64_t theBox = ii->second.id();
    unsigned theRank = NaluEnv::self().parallel_rank();
    const unsigned pt_proc = ii->first.proc();
    const unsigned box_proc = ii->second.proc();
    if ( (box_proc == theRank) && (pt_proc != theRank) ) {

      // Send box to pt proc

      // find the element
      stk::mesh::Entity theElemMeshObj = bulk_data.get_entity(stk::topology::ELEMENT_RANK, theBox);
      if ( !(bulk_data.is_valid(theElemMeshObj)) )
        throw std::runtime_error("no valid entry for element");

      // new element to ghost counter
      realm_.contactManager_->needToGhostCount_++;

      // deal with elements to push back to be ghosted
      stk::mesh::EntityProc theElemPair(theElemMeshObj, pt_proc);
      realm_.contactManager_->elemsToGhost_.push_back(theElemPair);
    }
  }
}

//--------------------------------------------------------------------------
//-------- complete_search -------------------------------------------------
//--------------------------------------------------------------------------
void
ContactInfo::complete_search()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  const int nDim = meta_data.spatial_dimension();

  // fields
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // extract master element; here note that we are only supporting hex or quad...
  const stk::topology & theElemTopo = (nDim == 2) ? stk::topology::QUAD_4_2D : stk::topology::HEX_8;
  meSCS_ = realm_.get_surface_master_element(theElemTopo);

  const int nodesPerElement = meSCS_->nodesPerElement_;

  std::vector<double> isoParCoords(nDim);
  std::vector<double> theElementCoords(nDim*nodesPerElement);
  std::vector<double> theHaloCoords(nDim);

  // now proceed with the standard search
  std::vector<std::pair<boundingPoint::second_type, boundingElementBox::second_type> >::const_iterator ii;
  for( ii=searchKeyPair_.begin(); ii!=searchKeyPair_.end(); ++ii ) {

    const uint64_t thePt = ii->first.id();
    const uint64_t theBox = ii->second.id();
    const unsigned theRank = NaluEnv::self().parallel_rank();
    const unsigned pt_proc = ii->first.proc();

    // check if I own the point...
    if ( theRank == pt_proc ) {

      // yes, I own the point... However, what about the element? Who owns that?

      // proceed as required; all elements should have already been ghosted via the coarse search
      stk::mesh::Entity elem = bulk_data.get_entity(stk::topology::ELEMENT_RANK, theBox);
      if ( !(bulk_data.is_valid(elem)) )
        throw std::runtime_error("no valid entry for element");

      int elemIsGhosted = bulk_data.bucket(elem).owned() ? 0 : 1;

      std::map<uint64_t, HaloInfo *>::iterator iterHalo;
      iterHalo=haloInfoMap_.find(thePt);
      if ( iterHalo == haloInfoMap_.end() )
        throw std::runtime_error("no valid entry for haloInfMap");

      HaloInfo *theHaloInfo = iterHalo->second;

      theHaloCoords = theHaloInfo->haloNodalCoords_;

      // now load the elemental nodal coords
      stk::mesh::Entity const * elem_node_rels = bulk_data.begin_nodes(elem);
      int num_nodes = bulk_data.num_nodes(elem);

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];
        const double * coords =  stk::mesh::field_data(*coordinates, node);
        for ( int j = 0; j < nDim; ++j ) {
          const int offSet = j*nodesPerElement +ni;
          theElementCoords[offSet] = coords[j];
        }
      }

      const double nearestDistance = meSCS_->isInElement(&theElementCoords[0],
                                                         &(theHaloCoords[0]),
                                                         &(isoParCoords[0]));
      if ( nearestDistance < theHaloInfo->bestX_ ) {
        theHaloInfo->owningElement_ = elem;
        theHaloInfo->isoParCoords_ = isoParCoords;
        theHaloInfo->bestX_ = nearestDistance;
        theHaloInfo->elemIsGhosted_ = elemIsGhosted;
      }
    }
    else {
      // not this proc's issue
    }

  }

  // check to see that all halo coords have a home...
  std::vector<double> haloCoordCheck(nDim);
  const double tol = 1.0e-6;
  const double maxTol = 1.0+tol;
  std::map<uint64_t, HaloInfo *>::iterator iterHalo;
  for (iterHalo  = haloInfoMap_.begin();
       iterHalo != haloInfoMap_.end();
       ++iterHalo) {

    // get halo info object and element
    HaloInfo * infoObject = (*iterHalo).second;
    stk::mesh::Entity elem = infoObject->owningElement_;
    if ( infoObject->bestX_ > maxTol || !(bulk_data.is_valid(elem)) ) {

      NaluEnv::self().naluOutputP0() << "Sorry, halo node for face node " << bulk_data.identifier(infoObject->faceNode_)
                      << " Does not have an ideal bestX; consider clipping "
                      << infoObject->bestX_ << std::endl;
      if ( !(bulk_data.is_valid(elem)) )
        NaluEnv::self().naluOutputP0() << "In fact, the owning master element is null" << std::endl;
      else
        NaluEnv::self().naluOutputP0() << " The best element is "
                        << bulk_data.identifier(infoObject->owningElement_)
                        << " consider clipping "<< std::endl;

      // clip to isoPar min/max... Still let dxj be what it is...
      if ( infoObject->bestX_ > maxTol && clipIsoParametricCoords_ ) {
        NaluEnv::self().naluOutputP0()
          << "Will clip the isoParametricCoords for node id: " << bulk_data.identifier(infoObject->faceNode_) << std::endl;
        const double minTol = -1.0-tol;
        for ( int j = 0; j < nDim; ++j ) {
          if ( infoObject->isoParCoords_[j] > maxTol )
            infoObject->isoParCoords_[j] = maxTol;
          if ( infoObject->isoParCoords_[j] < minTol)
            infoObject->isoParCoords_[j] = minTol;
        }
      }
    }
    else {
      // check to see if the isoparametric coords provides the coords of the halo node..
      theHaloCoords = infoObject->haloNodalCoords_;
      isoParCoords = infoObject->isoParCoords_;

      // now load the elemental nodal coords
      stk::mesh::Entity const * elem_node_rels = bulk_data.begin_nodes(elem);
      int num_nodes = bulk_data.num_nodes(elem);

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];
        const double * coords =  stk::mesh::field_data(*coordinates, node);
        for ( int j = 0; j < nDim; ++j ) {
          const int offSet = j*nodesPerElement +ni;
          theElementCoords[offSet] = coords[j];
        }
      }

      meSCS_->interpolatePoint(nDim,
                               &(isoParCoords[0]),
                               &theElementCoords[0],
                               &(haloCoordCheck[0]));
      for (int j = 0; j < nDim; ++j )
        infoObject->checkhaloNodalCoords_[j] = haloCoordCheck[j];
    }
  }
}

//--------------------------------------------------------------------------
//-------- find_possible_elements ------------------------------------------
//--------------------------------------------------------------------------
void
ContactInfo::find_possible_elements()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

  // fields
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // these two hacks are a consequence of contact not really
  // supporting 2d naturally
  Point minCorner, maxCorner;

  stk::mesh::Selector s_locally_owned = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(contactSearchBlock_);

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
      stk::mesh::Entity const* elem_node_rels = bulk_data.begin_nodes(elem);
      const int num_nodes = bulk_data.num_nodes(elem);

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

      // now loop over max min
      double rMax = 0.0;
      double rMin = 0.0;
      for ( int j = 0; j < nDim; ++j ) {
        rMax += minCorner[j]*minCorner[j];
        rMin += maxCorner[j]*maxCorner[j];
      }
      rMax = std::sqrt(rMax);
      rMin = std::sqrt(rMin);

      // check if this fits the bill...
      if ( (rMin <= maxSearchRadius_ && rMin >= minSearchRadius_ ) ||
           (rMax <= maxSearchRadius_ && rMax >= minSearchRadius_ )  ) {

        // setup ident
        stk::search::IdentProc<uint64_t,int> theIdent(bulk_data.identifier(elem), NaluEnv::self().parallel_rank());

        // expand the box
        for ( int i = 0; i < nDim; ++i ) {
          const double theMin = minCorner[i];
          const double theMax = maxCorner[i];
          const double increment = expandBoxPercentage_*(theMax - theMin);
          minCorner[i]   -= increment;
          maxCorner[i] += increment;
        }

        // create the bounding point box and push back
        boundingElementBox theBox(Box(minCorner,maxCorner), theIdent);
        boundingElementBoxVec_.push_back(theBox);
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- set_best_x ------------------------------------------------------
//--------------------------------------------------------------------------
void
ContactInfo::set_best_x()
{
  // show halo info
  std::map<uint64_t, HaloInfo *>::iterator iterHalo;
  for (iterHalo  = haloInfoMap_.begin();
       iterHalo != haloInfoMap_.end();
       ++iterHalo) {

    HaloInfo * infoObject = (*iterHalo).second;
    infoObject->bestX_ = 1.0e16;
  }
}

//--------------------------------------------------------------------------
//-------- dump_diagnosis --------------------------------------------------
//--------------------------------------------------------------------------
void
ContactInfo::dump_diagnosis()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

  std::map<uint64_t,HaloInfo *>::iterator iterHalo;
  for (iterHalo  = haloInfoMap_.begin();
       iterHalo != haloInfoMap_.end();
       ++iterHalo) {

    HaloInfo * infoObject = (*iterHalo).second;
    const size_t theId = bulk_data.identifier(infoObject->faceNode_);
    NaluEnv::self().naluOutputP0() << " " << std::endl;
    NaluEnv::self().naluOutputP0() << "The node id " << theId << " Has the following information: "
                    << infoObject->haloEdgeDs_ << std::endl;
    double norm = 0.0;
    for ( int j = 0; j < nDim; ++j ) {
      NaluEnv::self().naluOutputP0() << "direction: " << j << " "
                      << infoObject->haloEdgeAreaVec_[j] << " "
                      << infoObject->nodalCoords_[j] << " "
                      << infoObject->haloNodalCoords_[j] << " "
                      << infoObject->checkhaloNodalCoords_[j] << std::endl;
      const double diff = infoObject->haloNodalCoords_[j] - infoObject->checkhaloNodalCoords_[j];
      norm += diff*diff;
    }
    norm = std::sqrt(norm);
    if ( norm > 1.0e-8 )
      NaluEnv::self().naluOutputP0() << "WARNING: diference between halo coords and check halo coords is large " << norm << std::endl;

    if ( !(bulk_data.is_valid(infoObject->owningElement_)) ) {
      NaluEnv::self().naluOutputP0() << "There was NO ELEMENT FOUND!! Think about expanding the box or search band" << std::endl;
    }
    else {
      NaluEnv::self().naluOutputP0() << "Finally, the found element is " << bulk_data.identifier(infoObject->owningElement_)
                      << " with bestX as: " << infoObject->bestX_ << " " << std::endl;
    }

    const int isGhosted = infoObject->elemIsGhosted_;
    if ( isGhosted > 0 ) {
      NaluEnv::self().naluOutputP0() << "The node id " << theId << " found a need to ghost the element "
                      << bulk_data.identifier(infoObject->owningElement_) << std::endl;
    }
  }

  for (iterHalo  = haloInfoMap_.begin();
       iterHalo != haloInfoMap_.end();
       ++iterHalo) {

    HaloInfo * infoObject = (*iterHalo).second;
    const int theId = bulk_data.identifier(infoObject->faceNode_);
    const int isGhosted = infoObject->elemIsGhosted_;
    if ( isGhosted > 0 ) {
      NaluEnv::self().naluOutputP0() << "The node id " << theId << " found a need to ghost the element "
                      << bulk_data.identifier(infoObject->owningElement_) << std::endl;
    }
  }
}

} // namespace nalu
} // namespace sierra
