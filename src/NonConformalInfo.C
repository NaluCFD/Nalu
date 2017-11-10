/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <NonConformalInfo.h>
#include <NonConformalManager.h>
#include <DgInfo.h>
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

// stk_search
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>

// stk_topo
#include <stk_topology/topology.hpp>

// vector, pair and find
#include <vector>
#include <utility>
#include <algorithm>

namespace sierra{
namespace nalu{

class DgInfo;

// compare operator 
struct compareGaussPoint {
  compareGaussPoint()  {}
  bool operator () (const std::pair<theKey, theKey> &p, const uint64_t i) {
    return (p.first.id() < i);
  }
  bool operator () (const uint64_t i, const std::pair<theKey, theKey> &p) {
    return (i < p.first.id());
  }
};

// Sortintlowhigh
struct sortIntLowHigh {
  sortIntLowHigh() {}
  bool operator () (const std::pair<theKey, theKey> &p1, const std::pair<theKey, theKey> &p2) {
    return (p1.first.id() < p2.first.id());
  }
};

  
//==========================================================================
// Class Definition
//==========================================================================
// Acon_NonConformalInfo - contains virtual edge data
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
NonConformalInfo::NonConformalInfo(
   Realm &realm,
   const stk::mesh::PartVector currentPartVec,
   const stk::mesh::PartVector opposingPartVec,
   const double expandBoxPercentage,
   const std::string &searchMethodName,
   const bool clipIsoParametricCoords,
   const double searchTolerance,
   const std::string debugName)
  : realm_(realm ),
    name_(debugName),
    currentPartVec_(currentPartVec),
    opposingPartVec_(opposingPartVec),
    expandBoxPercentage_(expandBoxPercentage),
    searchMethod_(stk::search::BOOST_RTREE),
    clipIsoParametricCoords_(clipIsoParametricCoords),
    searchTolerance_(searchTolerance),
    meshMotion_(realm_.has_mesh_motion()),
    canReuse_(false)
{
  // determine search method for this pair
  if ( searchMethodName == "boost_rtree" )
    searchMethod_ = stk::search::BOOST_RTREE;
  else if ( searchMethodName == "stk_kdtree" )
    searchMethod_ = stk::search::KDTREE;
  else
    NaluEnv::self().naluOutputP0() << "NonConformalInfo::search method not declared; will use boost_rtree" << std::endl;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
NonConformalInfo::~NonConformalInfo()
{
  delete_dgInfo();
}

//--------------------------------------------------------------------------
//-------- delete_dgInfo ---------------------------------------------------
//--------------------------------------------------------------------------
void
NonConformalInfo::delete_dgInfo()
{
  std::vector<std::vector<DgInfo*> >::iterator ii;
  for( ii=dgInfoVec_.begin(); ii!=dgInfoVec_.end(); ++ii ) {
    std::vector<DgInfo *> &faceDgInfoVec = (*ii);
    for ( size_t k = 0; k < faceDgInfoVec.size(); ++k )
      delete faceDgInfoVec[k];
  }
  dgInfoVec_.clear();
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
NonConformalInfo::initialize()
{

  // clear some of the search info
  boundingPointVec_.clear();
  boundingElementBoxVec_.clear();
  searchKeyPair_.clear();

  // delete info only if adaptivity is active
  if ( realm_.mesh_changed() ) {
    delete_dgInfo();
  }
  
  // construct if the size is zero; reset always
  if ( dgInfoVec_.size() == 0 )
    construct_dgInfo();
  reset_dgInfo();
  
  // construct the points and boxes required for the search
  construct_bounding_points();
  construct_bounding_boxes();

  // ghosting
  determine_elems_to_ghost();
}

//--------------------------------------------------------------------------
//-------- construct_dgInfo ------------------------------------------------
//--------------------------------------------------------------------------
void
NonConformalInfo::construct_dgInfo()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

  // define vector of parent topos; should always be UNITY in size
  std::vector<stk::topology> parentTopo;
  
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(currentPartVec_);
  
  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
  
  // need to keep track of some sort of local id for each gauss point...
  uint64_t localGaussPointId = 0;
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    
    stk::mesh::Bucket & b = **ib;
    
    const stk::mesh::Bucket::size_type length   = b.size();
    
    // extract connected element topology
    b.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
    ThrowAssert ( parentTopo.size() == 1 );
    stk::topology currentElemTopo = parentTopo[0];

    // volume and surface master element
    MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(currentElemTopo);
    MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());
    
    // master element-specific values
    const int numScsBip = meFC->numIntPoints_;

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get face, global and local id
      stk::mesh::Entity face = b[k];
      uint64_t globalFaceId = bulk_data.identifier(face);
      
      // extract the connected element to this exposed face; should be single in size!
      const stk::mesh::Entity* face_elem_rels = bulk_data.begin_elements(face);
      ThrowAssert( bulk_data.num_elements(face) == 1 );

      // get element; its face ordinal number
      stk::mesh::Entity element = face_elem_rels[0];
      const stk::mesh::ConnectivityOrdinal* face_elem_ords = bulk_data.begin_element_ordinals(face);
      const int currentFaceOrdinal = face_elem_ords[0];
      
      std::vector<DgInfo *> faceDgInfoVec(numScsBip);
      for ( int ip = 0; ip < numScsBip; ++ip ) { 
        DgInfo *dgInfo = new DgInfo(NaluEnv::self().parallel_rank(), globalFaceId, localGaussPointId++, ip, 
                                    face, element, currentFaceOrdinal, meFC, meSCS, currentElemTopo, nDim); 
        faceDgInfoVec[ip] = dgInfo;
      }
      
      // push them all back
      dgInfoVec_.push_back(faceDgInfoVec);
    }
  }
}

//--------------------------------------------------------------------------
//-------- reset_dgInfo ----------------------------------------------------
//--------------------------------------------------------------------------
void
NonConformalInfo::reset_dgInfo()
{
  std::vector<std::vector<DgInfo*> >::iterator ii;
  for( ii=dgInfoVec_.begin(); ii!=dgInfoVec_.end(); ++ii ) {
    std::vector<DgInfo *> &theVec = (*ii);    
    for ( size_t k = 0; k < theVec.size(); ++k ) {
      DgInfo *dgInfo = theVec[k];
      if ( !canReuse_ ) {
        dgInfo->allOpposingElementIdsOld_.clear();
        dgInfo->allOpposingElementIdsOld_ = dgInfo->allOpposingElementIds_;
      }
      // always reset bestX and opposing faceIDs for the upcoming search
      dgInfo->bestX_ = dgInfo->bestXRef_;
      dgInfo->allOpposingElementIds_.clear();
    }
  }
}
  
//--------------------------------------------------------------------------
//-------- construct_bounding_points ---------------------------------------
//--------------------------------------------------------------------------
void
NonConformalInfo::construct_bounding_points()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

  // dial in Gauss-labato points
  const bool useShifted = realm_.has_nc_gauss_labatto_quadrature();

  // hold the point location for integration points
  Point currentIpCoords;

  // nodal fields to gather
  std::vector<double> ws_face_coordinates;
  // master element
  std::vector<double> ws_face_shape_function;

  // fields
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  
  std::vector<std::vector<DgInfo*> >::iterator ii;
  for( ii=dgInfoVec_.begin(); ii!=dgInfoVec_.end(); ++ii ) {
    std::vector<DgInfo *> &theVec = (*ii);

    //=======================================================
    // all ips on this face use a common face master element 
    //            gather common operations once 
    //=======================================================
    const DgInfo *firstDgInfo = theVec[0];
    MasterElement *meFC = firstDgInfo->meFCCurrent_;
    
    // master element-specific values
    const int numScsBip = meFC->numIntPoints_;
    const int nodesPerFace = meFC->nodesPerElement_;
   
    // algorithm related; face
    ws_face_coordinates.resize(nodesPerFace*nDim);  
    ws_face_shape_function.resize(numScsBip*nodesPerFace);
    
    // pointers
    double *p_face_coordinates = &ws_face_coordinates[0];
    double *p_face_shape_function = &ws_face_shape_function[0];
    
    // populate shape function
    if ( useShifted )
      meFC->shifted_shape_fcn(&p_face_shape_function[0]);
    else
      meFC->shape_fcn(&p_face_shape_function[0]);
  
    // gather nodal data off of face
    stk::mesh::Entity const * face_node_rels = bulk_data.begin_nodes(firstDgInfo->currentFace_);
    const int num_face_nodes = bulk_data.num_nodes(firstDgInfo->currentFace_);
    
    // sanity check on num nodes
    ThrowAssert( num_face_nodes == nodesPerFace );
    for ( int ni = 0; ni < num_face_nodes; ++ni ) {
      stk::mesh::Entity node = face_node_rels[ni];
      double * coords = stk::mesh::field_data(*coordinates, node);
      for ( int j =0; j < nDim; ++j ) {
        p_face_coordinates[ni*nDim+j] = coords[j];
      }
    }
    
    // now loop over all ips on this face
    for ( size_t k = 0; k < theVec.size(); ++k ) {
      
      DgInfo *dgInfo = theVec[k];

      // local and current ip
      const uint64_t localIp  = dgInfo->localGaussPointId_; 
      const int currentFaceIp = dgInfo->currentGaussPointId_;

      // compute coordinates
      for ( int j = 0; j < nDim; ++j )
        currentIpCoords[j] = 0.0;
        
      // interpolate to gauss point
      for ( int ic = 0; ic < nodesPerFace; ++ic ) {
        const double r = p_face_shape_function[currentFaceIp*nodesPerFace+ic];
        for ( int j = 0; j < nDim; ++j ) {
          currentIpCoords[j] += r*p_face_coordinates[ic*nDim+j];
        }
      }

      // extract isoparametric coords on current face from meFC
      const double *intgLoc = useShifted ? &meFC->intgLocShift_[0] : &meFC->intgLoc_[0];
      
      // copy these coordinates
      for ( int j = 0; j < nDim; ++j ) {
        dgInfo->currentGaussPointCoords_[j] = currentIpCoords[j];
      }
      
      // save face iso-parametric coordinates; extract conversion factor from CVFEM to isInElement
      const double conversionFac = meFC->scaleToStandardIsoFac_;
      for ( int j = 0; j < nDim-1; ++j ) {
        dgInfo->currentIsoParCoords_[j] = conversionFac*intgLoc[currentFaceIp*(nDim-1)+j]; 
      }
      
      // setup ident for this point; use local integration point id
      stk::search::IdentProc<uint64_t,int> theIdent(localIp, NaluEnv::self().parallel_rank());
      
      // create the bounding point and push back
      boundingPoint thePt(currentIpCoords, theIdent);
      boundingPointVec_.push_back(thePt);
    } 
  }
}

//--------------------------------------------------------------------------
//-------- determine_elems_to_ghost ----------------------------------------
//--------------------------------------------------------------------------
void
NonConformalInfo::determine_elems_to_ghost()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // perform the coarse search
  stk::search::coarse_search(boundingPointVec_, boundingElementBoxVec_, searchMethod_, NaluEnv::self().parallel_comm(), searchKeyPair_);

  // sort based on local gauss point
  std::sort (searchKeyPair_.begin(), searchKeyPair_.end(), sortIntLowHigh());

  std::vector<std::pair<theKey, theKey> >::const_iterator ii;
  for( ii=searchKeyPair_.begin(); ii!=searchKeyPair_.end(); ++ii ) {

    const uint64_t theBox = ii->second.id();
    unsigned theRank = NaluEnv::self().parallel_rank();
    const unsigned pt_proc = ii->first.proc();
    const unsigned box_proc = ii->second.proc();
    if ( (box_proc == theRank) && (pt_proc != theRank) ) {

      // Send box to pt proc

      // find the face element
      stk::mesh::Entity element = bulk_data.get_entity(stk::topology::ELEMENT_RANK, theBox);
      if ( !(bulk_data.is_valid(element)) )
        throw std::runtime_error("no valid entry for element in determine_elems_to_ghost");

      // deal with elements to push back to be ghosted; downward relations come for the ride...
      stk::mesh::EntityProc theElemPair(element, pt_proc);
      realm_.nonConformalManager_->elemsToGhost_.push_back(theElemPair);
    }
  }
}

//--------------------------------------------------------------------------
//-------- complete_search -------------------------------------------------
//--------------------------------------------------------------------------
void
NonConformalInfo::complete_search()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  const int nDim = meta_data.spatial_dimension();

  // fields
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  std::vector<double> currentGaussPointCoords(nDim);
  std::vector<double> opposingIsoParCoords(nDim);

  // invert the process... Loop over dgInfoVec_ and query searchKeyPair_ for this information
  std::vector<DgInfo *> problemDgInfoVec;
  std::vector<std::vector<DgInfo*> >::iterator ii;
  for( ii=dgInfoVec_.begin(); ii!=dgInfoVec_.end(); ++ii ) {
    std::vector<DgInfo *> &theVec = (*ii);
    for ( size_t k = 0; k < theVec.size(); ++k ) {
      
      DgInfo *dgInfo = theVec[k];
      const uint64_t localGaussPointId  = dgInfo->localGaussPointId_; 

      std::pair <std::vector<std::pair<theKey, theKey> >::const_iterator, std::vector<std::pair<theKey, theKey> >::const_iterator > 
        p2 = std::equal_range(searchKeyPair_.begin(), searchKeyPair_.end(), localGaussPointId, compareGaussPoint());

      if ( p2.first == p2.second ) {
        problemDgInfoVec.push_back(dgInfo);        
      }
      else {
        for (std::vector<std::pair<theKey, theKey> >::const_iterator jj = p2.first; jj != p2.second; ++jj ) {
          
          const uint64_t theBox = jj->second.id();
          const unsigned theRank = NaluEnv::self().parallel_rank();
          const unsigned pt_proc = jj->first.proc();
          
          // check if I own the point...
          if ( theRank == pt_proc ) {
            
            // yes, I own the point... However, what about the face element? Who owns that?

            // proceed as required; all elements should have already been ghosted via the coarse search
            stk::mesh::Entity opposingElement = bulk_data.get_entity(stk::topology::ELEMENT_RANK, theBox);
            if ( !(bulk_data.is_valid(opposingElement)) )
              throw std::runtime_error("no valid entry for element");

            int opposingElementIsGhosted = bulk_data.bucket(opposingElement).owned() ? 0 : 1;
            
            // extract the gauss point coordinates
            currentGaussPointCoords = dgInfo->currentGaussPointCoords_;
            
            // now load the elemental nodal coords
            stk::mesh::Entity const * elem_node_rels = bulk_data.begin_nodes(opposingElement);
            int num_nodes = bulk_data.num_nodes(opposingElement);
            
            std::vector<double> theElementCoords(nDim*num_nodes);
            
            for ( int ni = 0; ni < num_nodes; ++ni ) {
              stk::mesh::Entity node = elem_node_rels[ni];
              const double * coords =  stk::mesh::field_data(*coordinates, node);
              for ( int j = 0; j < nDim; ++j ) {
                const int offSet = j*num_nodes +ni;
                theElementCoords[offSet] = coords[j];
              }
            }
                        
            // extract the opposing element topo and associated master element
            const stk::topology theOpposingElementTopo = bulk_data.bucket(opposingElement).topology();
            MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(theOpposingElementTopo);
            
            // possible reuse            
            dgInfo->allOpposingElementIds_.push_back(theBox);
            
            // find distance between true current gauss point coords (the point) and the candidate bounding box
            const double nearestDistance = meSCS->isInElement(&theElementCoords[0],
                                                              &(currentGaussPointCoords[0]),
                                                              &(opposingIsoParCoords[0]));
            
            // check is this is the best candidate
            if ( nearestDistance < dgInfo->bestX_ ) {
              // save off all required opposing information
              dgInfo->opposingElement_ = opposingElement;
              dgInfo->meSCSOpposing_ = meSCS;
              dgInfo->opposingElementTopo_ = theOpposingElementTopo;
              dgInfo->opposingIsoParCoords_ = opposingIsoParCoords;
              dgInfo->bestX_ = nearestDistance;
              dgInfo->opposingElementIsGhosted_ = opposingElementIsGhosted;
            }
          }
          else {
            // not this proc's issue
          }
        }
      }
    }
  }
  
  // check for problems... will want to be more pro-active in the near future, e.g., expand and search...
  if ( problemDgInfoVec.size() > 0 ) {
    NaluEnv::self().naluOutputP0() << "NonConformalInfo::complete_search issue with " << name_ 
                                   << " Size of issue is " << problemDgInfoVec.size() << std::endl; 
    NaluEnv::self().naluOutputP0() << "Problem ips are as follows: " << std::endl; 
    for ( size_t k = 0; k < problemDgInfoVec.size(); ++k ) {
      problemDgInfoVec[k]->dump_info(); 
    }
    NaluEnv::self().naluOutputP0() << std::endl;
    throw std::runtime_error("Try to adjust the search tolerance and re-submit...");
  }

  // check for reuse (debug now)
  size_t numberOfFacesMissing = 0;
  for( size_t iv = 0; iv < dgInfoVec_.size(); ++iv ) {
    std::vector<DgInfo *> &theVec = dgInfoVec_[iv];
    for ( size_t k = 0; k < theVec.size(); ++k ) {
      
      // extract the info object; new and old
      DgInfo *dgInfo = theVec[k];
      
      // extract the bestX opposing face id
      const size_t bestOpposingId = bulk_data.identifier(dgInfo->opposingElement_);
      
      // is the required active stencil opposing id within the vector of face 
      // ids returned in the search? (no need to sort given the size)
      auto itF = std::find(dgInfo->allOpposingElementIdsOld_.begin(), 
                           dgInfo->allOpposingElementIdsOld_.end(), bestOpposingId);
      
      // increment missing faces if NOT found
      if (itF == dgInfo->allOpposingElementIdsOld_.end()) {
        numberOfFacesMissing++;
      }
    }
  }
  
  // global sum
  size_t g_numberOfFacesMissing;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &numberOfFacesMissing, &g_numberOfFacesMissing, 1);
  if ( g_numberOfFacesMissing > 0 ) {
    NaluEnv::self().naluOutputP0() << "Ghosted search entries ARE NOT sufficient for re-use " << std::endl;
    canReuse_ = false;
  }
  else {
    NaluEnv::self().naluOutputP0() << "Ghosted search entries ARE sufficient for re-use " << std::endl;
    canReuse_ = false; // Set the false until ready for primetime
  }
}
  
//--------------------------------------------------------------------------
//-------- construct_bounding_boxes ----------------------------------------
//--------------------------------------------------------------------------
void
NonConformalInfo::construct_bounding_boxes()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

  // fields
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // points
  Point minCorner, maxCorner;

  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(opposingPartVec_);

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );

  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {

    stk::mesh::Bucket & b = **ib;

    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get face
      stk::mesh::Entity face = b[k];

      // extract the connected element to this exposed face; should be single in size!
      const stk::mesh::Entity* face_elem_rels = bulk_data.begin_elements(face);
      ThrowAssert( bulk_data.num_elements(face) == 1 );
      
      // get element; its face ordinal number and populate face_node_ordinals
      stk::mesh::Entity element = face_elem_rels[0];

      // initialize max and min
      for (int j = 0; j < nDim; ++j ) {
        minCorner[j] = +1.0e16;
        maxCorner[j] = -1.0e16;
      }

      // extract elem_node_relations
      stk::mesh::Entity const* elem_node_rels = bulk_data.begin_nodes(element);
      const int num_nodes = bulk_data.num_nodes(element);

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
      stk::search::IdentProc<uint64_t,int> theIdent(bulk_data.identifier(element), NaluEnv::self().parallel_rank());

      // expand the box by both % and search tolerance
      for ( int i = 0; i < nDim; ++i ) {
        const double theMin = minCorner[i];
        const double theMax = maxCorner[i];
        const double increment = expandBoxPercentage_*(theMax - theMin) + searchTolerance_;
        minCorner[i] -= increment;
        maxCorner[i] += increment;
      }

      // create the bounding point box and push back
      boundingElementBox theBox(Box(minCorner,maxCorner), theIdent);
      boundingElementBoxVec_.push_back(theBox);
    }
  }
}  

//--------------------------------------------------------------------------
//-------- provide_diagnosis -----------------------------------------------
//--------------------------------------------------------------------------
void
NonConformalInfo::provide_diagnosis()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  const int nDim = meta_data.spatial_dimension();

  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  std::vector<double> currentGaussPointCoords(nDim);
  std::vector<double> opposingGaussPointCoords(nDim);
 
  std::vector<double> currentIsoParCoords(nDim);
  std::vector<double> opposingIsoParCoords(nDim);

  NaluEnv::self().naluOutput() << std::endl;
  NaluEnv::self().naluOutput() << "Non Conformal Alg review for surface: " << name_ << std::endl;
  NaluEnv::self().naluOutput() << "===================================== " << std::endl;
  std::vector<std::vector<DgInfo*> >::iterator ii;
  for( ii=dgInfoVec_.begin(); ii!=dgInfoVec_.end(); ++ii ) {
    std::vector<DgInfo *> &theVec = (*ii);
    for ( size_t k = 0; k < theVec.size(); ++k ) {
      DgInfo *dgInfo = theVec[k];

      // first, dump info
      dgInfo->dump_info();

      // now proceed to detailed face/element current/opposing checks
      const uint64_t localGaussPointId  = dgInfo->localGaussPointId_; 
      const uint64_t currentGaussPointId  = dgInfo->currentGaussPointId_; 

      // extract current face
      stk::mesh::Entity currentFace = dgInfo->currentFace_;

      // extract the gauss point isopar/geometric coordinates for current
      currentGaussPointCoords = dgInfo->currentGaussPointCoords_;
      currentIsoParCoords = dgInfo->currentIsoParCoords_;

      // extract the master element for current; with npe
      MasterElement *meFCCurrent = dgInfo->meFCCurrent_;      
      const int currentNodesPerFace = meFCCurrent->nodesPerElement_;

      // face:node relations
      stk::mesh::Entity const * current_face_node_rels = bulk_data.begin_nodes(currentFace);
      int current_face_num_nodes = bulk_data.num_nodes(currentFace); 

      // gather nodal coordinates
      std::vector <double > currentFaceNodalCoords(nDim*currentNodesPerFace);
      for ( int ni = 0; ni < current_face_num_nodes; ++ni ) {
        stk::mesh::Entity node = current_face_node_rels[ni];
        const double * coords = stk::mesh::field_data(*coordinates, node);
        for ( int j=0; j < nDim; ++j ) {
          currentFaceNodalCoords[j*currentNodesPerFace+ni] = coords[j];
        }
      }

      // interpolate to current GP
      std::vector<double> checkCurrentFaceGaussPointCoords(nDim);
      meFCCurrent->interpolatePoint(
        nDim,
        &currentIsoParCoords[0],
        &currentFaceNodalCoords[0],
        &checkCurrentFaceGaussPointCoords[0]);
      
      // extract current element
      stk::mesh::Entity currentElement = dgInfo->currentElement_;
      
      // best X
      const double bX = dgInfo->bestX_ ;
      
      // best opposing element
      stk::mesh::Entity theBestElement = dgInfo->opposingElement_;
      
      // extract the gauss point isopar coordiantes for opposing
      opposingIsoParCoords = dgInfo->opposingIsoParCoords_;

      // extract the master element for opposing; with npe
      MasterElement *meSCSOpposing = dgInfo->meSCSOpposing_;      
      const int opposingNodesPerElement = meSCSOpposing->nodesPerElement_;

      // face:node relations
      stk::mesh::Entity const * opposing_element_node_rels = bulk_data.begin_nodes(theBestElement);
      int opposing_element_num_nodes = bulk_data.num_nodes(theBestElement); 

      // gather nodal coordinates
      std::vector <double > opposingElementNodalCoords(nDim*opposingNodesPerElement);
      for ( int ni = 0; ni < opposing_element_num_nodes; ++ni ) {
        stk::mesh::Entity node = opposing_element_node_rels[ni];
        const double * coords = stk::mesh::field_data(*coordinates, node);
        for ( int j=0; j < nDim; ++j ) {
          opposingElementNodalCoords[j*opposingNodesPerElement+ni] = coords[j];
        }
      }
      
      // interpolate to opposing GP
      std::vector<double> checkOpposingElementGaussPointCoords(nDim);
      meSCSOpposing->interpolatePoint(
        nDim,
        &opposingIsoParCoords[0],
        &opposingElementNodalCoords[0],
        &(checkOpposingElementGaussPointCoords[0]));

      // global id for opposing element
      const uint64_t opElementId = bulk_data.identifier(dgInfo->opposingElement_);

      // compute a norm between the curent nd opposing coordinate checks
      double distanceNorm = 0.0;
      for ( int j = 0; j < nDim; ++j ) 
        distanceNorm += std::pow(checkCurrentFaceGaussPointCoords[j] -checkOpposingElementGaussPointCoords[j], 2);
      distanceNorm = std::sqrt(distanceNorm);

      // provide output...
      NaluEnv::self().naluOutput() << "Gauss Point Lid: " << localGaussPointId << " Review " << std::endl;
      NaluEnv::self().naluOutput() << "  encapsulated by Gid: (";  
      for ( int ni = 0; ni < current_face_num_nodes; ++ni ) {
        stk::mesh::Entity node = current_face_node_rels[ni];
        NaluEnv::self().naluOutput() << bulk_data.identifier(node) << " ";
      }
      NaluEnv::self().naluOutput() << ")" << std::endl;
      NaluEnv::self().naluOutput() << "Current Gauss Point id: " << currentGaussPointId 
                                   << " (nearest node) " << bulk_data.identifier(current_face_node_rels[currentGaussPointId])
                                   << std::endl;
      
      NaluEnv::self().naluOutput() << "  Current element Gid: " << bulk_data.identifier(currentElement) 
                                   << " (face ordinal: " << dgInfo->currentFaceOrdinal_ << ")" << std::endl;  
      
      NaluEnv::self().naluOutput() << "  has Gp coordinates: " ;
      for ( int i = 0; i < nDim; ++i )
        NaluEnv::self().naluOutput() << currentGaussPointCoords[i] << " ";
      NaluEnv::self().naluOutput() << std::endl;
      NaluEnv::self().naluOutput() << "  The best X is: " << bX << std::endl;
      NaluEnv::self().naluOutput() << "  Opposing element Gid: " << opElementId << std::endl;  
      NaluEnv::self().naluOutput() << "  encapsulated by Gid: (";  
      for ( int ni = 0; ni < opposing_element_num_nodes; ++ni ) {
        stk::mesh::Entity node = opposing_element_node_rels[ni];
        NaluEnv::self().naluOutput() << bulk_data.identifier(node) << " ";
      }
      NaluEnv::self().naluOutput() << ")" << std::endl;
      NaluEnv::self().naluOutput() << "  INTERNAL CHECK.... does current Gp and opposing found Gp match coordiantes? What error?" << std::endl;
      NaluEnv::self().naluOutput() << "  current and opposing Gp coordinates:        " << std::endl;
      for ( int i = 0; i < nDim; ++i )
        NaluEnv::self().naluOutput() << "      " << i << " " << checkCurrentFaceGaussPointCoords[i] << " " << checkOpposingElementGaussPointCoords[i] << std::endl;
      NaluEnv::self().naluOutput() << "  current and opposing Gp isoPar coordinates: " << std::endl;
      for ( int i = 0; i < nDim-1; ++i )
        NaluEnv::self().naluOutput() << "      " << i << " " << currentIsoParCoords[i] << " " << opposingIsoParCoords[i] << std::endl;
      NaluEnv::self().naluOutput() << std::endl;
      NaluEnv::self().naluOutput() << " in the end, the Error Distance Norm is: " << distanceNorm << std::endl;
      NaluEnv::self().naluOutput() << "-------------------------------------------------------------------" << std::endl;
    }
  }
}

//--------------------------------------------------------------------------
//-------- error_check -----------------------------------------------------
//--------------------------------------------------------------------------
size_t
NonConformalInfo::error_check()
{
  // check for coincident nodes via intersection of parts provided
  std::vector<stk::mesh::EntityId> coindidentNodesVec;

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  stk::mesh::Selector s_locally_owned_intersected = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(currentPartVec_) 
    &stk::mesh::selectUnion(opposingPartVec_);
  
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_locally_owned_intersected );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      coindidentNodesVec.push_back(bulk_data.identifier(b[k]));
    }
  } 
   
  // report the data if problem nodes were found
  if ( coindidentNodesVec.size() > 0 ) {
    NaluEnv::self().naluOutput() << std::endl;
    NaluEnv::self().naluOutput() << "Non Conformal Alg (P" << NaluEnv::self().parallel_rank() << ") error found on surface: " 
                                 << name_ << std::endl;
    NaluEnv::self().naluOutput() << "========================================= " << std::endl;
    for ( size_t k = 0; k < coindidentNodesVec.size(); ++k )
      NaluEnv::self().naluOutput() << "coincident nodeId found: " << coindidentNodesVec[k] << std::endl;
  }

  return coindidentNodesVec.size();
}

} // namespace nalu
} // namespace sierra
