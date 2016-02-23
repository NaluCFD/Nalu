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
#include <stk_util/environment/CPUTime.hpp>

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
   const stk::mesh::Part *currentPart,
   const stk::mesh::Part *opposingPart,
   const double expandBoxPercentage,
   const std::string &searchMethodName,
   const bool clipIsoParametricCoords,
   const double searchTolerance)
  : realm_(realm ),
    name_(currentPart->name()),
    currentPart_(currentPart),
    opposingPart_(opposingPart),
    expandBoxPercentage_(expandBoxPercentage),
    searchMethod_(stk::search::BOOST_RTREE),
    clipIsoParametricCoords_(clipIsoParametricCoords),
    searchTolerance_(searchTolerance),
    meshMotion_(realm_.has_mesh_motion())
{
  // determine search method for this pair
  if ( searchMethodName == "boost_rtree" )
    searchMethod_ = stk::search::BOOST_RTREE;
  else if ( searchMethodName == "stk_octree" )
    searchMethod_ = stk::search::OCTREE;
  else
    NaluEnv::self().naluOutputP0() << "NonConformalInfo::search method not declared; will use BOOST_RTREE" << std::endl;

}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
NonConformalInfo::~NonConformalInfo()
{
  // delete dgInfo info objects:
  std::vector<std::vector<DgInfo*> >::iterator ii;
  for( ii=dgInfoVec_.begin(); ii!=dgInfoVec_.end(); ++ii ) {
    std::vector<DgInfo *> &faceDgInfoVec = (*ii);
    for ( size_t k = 0; k < faceDgInfoVec.size(); ++k )
      delete faceDgInfoVec[k];
  }
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
NonConformalInfo::initialize()
{

  // clear some of the search info
  boundingPointVec_.clear();
  boundingFaceElementBoxVec_.clear();
  searchKeyPair_.clear();

  // delete dgInfoVec_
  std::vector<std::vector<DgInfo *> >::iterator ii;
  for( ii=dgInfoVec_.begin(); ii!=dgInfoVec_.end(); ++ii ) {
    std::vector<DgInfo *> &faceDgInfoVec = (*ii);
    for ( size_t k = 0; k < faceDgInfoVec.size(); ++k )
      delete faceDgInfoVec[k];
  }
  dgInfoVec_.clear();

  construct_dgInfo_state();

  find_possible_face_elements();

  determine_elems_to_ghost();

}

//--------------------------------------------------------------------------
//-------- construct_dgInfo_state --------------------------------------------
//--------------------------------------------------------------------------
void
NonConformalInfo::construct_dgInfo_state()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

  // dial in Gauss-labato points
  const bool useShifted = realm_.has_nc_gauss_labatto_quadrature();

  // hold the point location for gauss points
  Point currentGaussPointCoords;

  // nodal fields to gather
  std::vector<double> ws_face_coordinates;
  // master element
  std::vector<double> ws_face_shape_function;

  // define vector of parent topos; should always be UNITY in size
  std::vector<stk::topology> parentTopo;

  // fields
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::Selector(*currentPart_);

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
    MasterElement *meSCS = realm_.get_surface_master_element(currentElemTopo);
    MasterElement *meFC = realm_.get_surface_master_element(b.topology());

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

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get face, global and local id
      stk::mesh::Entity face = b[k];
      uint64_t globalFaceId = bulk_data.identifier(face);

      //======================================
      // gather nodal data off of face
      //======================================
      stk::mesh::Entity const * face_node_rels = bulk_data.begin_nodes(face);
      const int num_face_nodes = bulk_data.num_nodes(face);
      
      // sanity check on num nodes (low order, P=1, check)
      ThrowAssert( num_face_nodes == numScsBip ); ThrowAssert( num_face_nodes == nodesPerFace );
      for ( int ni = 0; ni < num_face_nodes; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        double * coords = stk::mesh::field_data(*coordinates, node);
        for ( int j =0; j < nDim; ++j ) {
          p_face_coordinates[ni*nDim+j] = coords[j];
        }
      }

      // extract the connected element to this exposed face; should be single in size!
      const stk::mesh::Entity* face_elem_rels = bulk_data.begin_elements(face);
      ThrowAssert( bulk_data.num_elements(face) == 1 );

      // get element; its face ordinal number
      stk::mesh::Entity element = face_elem_rels[0];
      const stk::mesh::ConnectivityOrdinal* face_elem_ords = bulk_data.begin_element_ordinals(face);
      const int currentFaceOrdinal = face_elem_ords[0];

      std::vector<DgInfo *> faceDgInfoVec(numScsBip);
      for ( int ip = 0; ip < numScsBip; ++ip ) {
        
        for ( int j = 0; j < nDim; ++j )
          currentGaussPointCoords[j] = 0.0;
        
        // interpolate to gauss point
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_face_shape_function[ip*nodesPerFace+ic];
          for ( int j = 0; j < nDim; ++j ) {
            currentGaussPointCoords[j] += r*p_face_coordinates[ic*nDim+j];
          }
        }
     
        // create data structure to hold this information; add currentIpNumber for later fast look-up
        DgInfo *dgInfo = new DgInfo(NaluEnv::self().parallel_rank(), globalFaceId, localGaussPointId, ip, 
                                    face, element, currentFaceOrdinal, meFC, meSCS, currentElemTopo, nDim);

        // extract isoparametric coords on current face from meFC
        const double *intgLoc = useShifted ? &meFC->intgLocShift_[0] : &meFC->intgLoc_[0];

        // copy these coordinates
        for ( int j = 0; j < nDim; ++j ) {
          dgInfo->currentGaussPointCoords_[j] = currentGaussPointCoords[j];
        }

        // save face iso-parametric coordinates; extract conversion factor from CVFEM to isInElement
        const double conversionFac = meFC->scaleToStandardIsoFac_;
        for ( int j = 0; j < nDim-1; ++j ) {
          dgInfo->currentIsoParCoords_[j] = conversionFac*intgLoc[ip*(nDim-1)+j]; 
        }

        // push back to local
        faceDgInfoVec[ip] = dgInfo;

        // setup ident for this point; use local gauss point id
        stk::search::IdentProc<uint64_t,int> theIdent(localGaussPointId++, NaluEnv::self().parallel_rank());

        // create the bounding point and push back
        boundingPoint thePt(currentGaussPointCoords, theIdent);
        boundingPointVec_.push_back(thePt);
      }
      
      // push them all back
      dgInfoVec_.push_back(faceDgInfoVec);
      
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
  stk::search::coarse_search(boundingPointVec_, boundingFaceElementBoxVec_, searchMethod_, NaluEnv::self().parallel_comm(), searchKeyPair_);

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
      stk::mesh::Entity face = bulk_data.get_entity(meta_data.side_rank(), theBox);
      if ( !(bulk_data.is_valid(face)) )
        throw std::runtime_error("no valid entry for face");
      
      // extract the connected element
      const stk::mesh::Entity* face_elem_rels = bulk_data.begin_elements(face);
      ThrowAssert( bulk_data.num_elements(face) == 1 );
      stk::mesh::Entity element = face_elem_rels[0];
          
      // new element to ghost counter
      realm_.nonConformalManager_->needToGhostCount_++;

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
        for (std::vector<std::pair<theKey, theKey> >::const_iterator ii = p2.first; ii != p2.second; ++ii ) {
          
          const uint64_t theBox = ii->second.id();
          const unsigned theRank = NaluEnv::self().parallel_rank();
          const unsigned pt_proc = ii->first.proc();
          
          // check if I own the point...
          if ( theRank == pt_proc ) {
            
            // yes, I own the point... However, what about the face element? Who owns that?

            // proceed as required; all elements should have already been ghosted via the coarse search
            stk::mesh::Entity opposingFace = bulk_data.get_entity(meta_data.side_rank(), theBox);
            if ( !(bulk_data.is_valid(opposingFace)) )
              throw std::runtime_error("no valid entry for face element");

            int opposingFaceIsGhosted = bulk_data.bucket(opposingFace).owned() ? 0 : 1;
            
            // extract the gauss point coordinates
            currentGaussPointCoords = dgInfo->currentGaussPointCoords_;
            
            // now load the face elemental nodal coords
            stk::mesh::Entity const * face_node_rels = bulk_data.begin_nodes(opposingFace);
            int num_nodes = bulk_data.num_nodes(opposingFace);
            
            std::vector<double> theElementCoords(nDim*num_nodes);
            
            for ( int ni = 0; ni < num_nodes; ++ni ) {
              stk::mesh::Entity node = face_node_rels[ni];
              const double * coords =  stk::mesh::field_data(*coordinates, node);
              for ( int j = 0; j < nDim; ++j ) {
                const int offSet = j*num_nodes +ni;
                theElementCoords[offSet] = coords[j];
              }
            }
            
            // extract the topo from this face element...
            const stk::topology theFaceTopo = bulk_data.bucket(opposingFace).topology();
            MasterElement *meFC = realm_.get_surface_master_element(theFaceTopo);
            
            // find distance between true current gauss point coords (the point) and the candidate bounding box
            const double nearestDistance = meFC->isInElement(&theElementCoords[0],
                                                             &(currentGaussPointCoords[0]),
                                                             &(opposingIsoParCoords[0]));
            if ( nearestDistance < dgInfo->bestX_ ) {
              // save the opposing face element and master element
              dgInfo->opposingFace_ = opposingFace;
              dgInfo->meFCOpposing_ = meFC;
              
              // extract the connected element to the opposing face
              const stk::mesh::Entity* face_elem_rels = bulk_data.begin_elements(opposingFace);
              ThrowAssert( bulk_data.num_elements(opposingFace) == 1 );
              stk::mesh::Entity opposingElement = face_elem_rels[0];
              dgInfo->opposingElement_ = opposingElement;

              // save off ordinal for opposing face
              const stk::mesh::ConnectivityOrdinal* face_elem_ords = bulk_data.begin_element_ordinals(opposingFace);
              dgInfo->opposingFaceOrdinal_ = face_elem_ords[0];
              
              // extract the opposing element topo and associated master element
              const stk::topology theOpposingElementTopo = bulk_data.bucket(opposingElement).topology();
              MasterElement *meSCS = realm_.get_surface_master_element(theOpposingElementTopo);
              dgInfo->meSCSOpposing_ = meSCS;
              dgInfo->opposingElementTopo_ = theOpposingElementTopo;
              dgInfo->opposingIsoParCoords_ = opposingIsoParCoords;
              dgInfo->bestX_ = nearestDistance;
              dgInfo->opposingFaceIsGhosted_ = opposingFaceIsGhosted;
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
    NaluEnv::self().naluOutputP0() << "NonConformalInfo::complete_search issue with " << currentPart_->name() << " " << opposingPart_->name() << " Size of issue is " << problemDgInfoVec.size() << std::endl; 
    NaluEnv::self().naluOutputP0() << "Problem ips are as follows: " << std::endl; 
    for ( size_t k = 0; k < problemDgInfoVec.size(); ++k ) {
      const uint64_t localGaussPointId  = problemDgInfoVec[k]->localGaussPointId_; 
      NaluEnv::self().naluOutputP0() << "local gauss point id with gass point coords " << localGaussPointId << " ";
      for ( int i = 0; i < nDim; ++i ) 
        NaluEnv::self().naluOutputP0() << " " << problemDgInfoVec[k]->currentGaussPointCoords_[i];
      NaluEnv::self().naluOutputP0() << std::endl;
    }
    NaluEnv::self().naluOutputP0() << std::endl;
    throw std::runtime_error("Try to adjust the search tolerance and re-submit...");
  }
}

//--------------------------------------------------------------------------
//-------- find_possible_face_elements -------------------------------------
//--------------------------------------------------------------------------
void
NonConformalInfo::find_possible_face_elements()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

  // fields
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // points
  Point minCorner, maxCorner;

  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::Selector(*opposingPart_);

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );

  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {

    stk::mesh::Bucket & b = **ib;

    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get face
      stk::mesh::Entity face = b[k];

      // initialize max and min
      for (int j = 0; j < nDim; ++j ) {
        minCorner[j] = +1.0e16;
        maxCorner[j] = -1.0e16;
      }

      // extract elem_node_relations
      stk::mesh::Entity const* face_node_rels = bulk_data.begin_nodes(face);
      const int num_nodes = bulk_data.num_nodes(face);

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];

        // pointers to real data
        const double * coords = stk::mesh::field_data(*coordinates, node );

        // check max/min
        for ( int j = 0; j < nDim; ++j ) {
          minCorner[j] = std::min(minCorner[j], coords[j]);
          maxCorner[j] = std::max(maxCorner[j], coords[j]);
        }
      }
      
      // setup ident
      stk::search::IdentProc<uint64_t,int> theIdent(bulk_data.identifier(face), NaluEnv::self().parallel_rank());

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
      boundingFaceElementBoxVec_.push_back(theBox);
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
      
      // best opposing face
      stk::mesh::Entity theBestFace = dgInfo->opposingFace_;
      
      // extract the gauss point isopar coordiantes for opposing
      opposingIsoParCoords = dgInfo->opposingIsoParCoords_;

      // extract the master element for opposing; with npe
      MasterElement *meFCOpposing = dgInfo->meFCOpposing_;      
      const int opposingNodesPerFace = meFCOpposing->nodesPerElement_;

      // face:node relations
      stk::mesh::Entity const * opposing_face_node_rels = bulk_data.begin_nodes(theBestFace);
      int opposing_face_num_nodes = bulk_data.num_nodes(theBestFace); 

      // gather nodal coordinates
      std::vector <double > opposingFaceNodalCoords(nDim*opposingNodesPerFace);
      for ( int ni = 0; ni < opposing_face_num_nodes; ++ni ) {
        stk::mesh::Entity node = opposing_face_node_rels[ni];
        const double * coords = stk::mesh::field_data(*coordinates, node);
        for ( int j=0; j < nDim; ++j ) {
          opposingFaceNodalCoords[j*opposingNodesPerFace+ni] = coords[j];
        }
      }
      
      // interpolate to opposing GP
      std::vector<double> checkOpposingFaceGaussPointCoords(nDim);
      meFCOpposing->interpolatePoint(
        nDim,
        &opposingIsoParCoords[0],
        &opposingFaceNodalCoords[0],
        &(checkOpposingFaceGaussPointCoords[0]));

      // global id for opposing element
      const uint64_t opElemId = bulk_data.identifier(dgInfo->opposingElement_);

      // compute a norm between the curent nd opposing coordinate checks
      double distanceNorm = 0.0;
      for ( int j = 0; j < nDim; ++j ) 
        distanceNorm += std::pow(checkCurrentFaceGaussPointCoords[j] -checkOpposingFaceGaussPointCoords[j], 2);
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
      NaluEnv::self().naluOutput() << "  Opposing element Gid: " << opElemId 
                                   << " (face ordinal: " << dgInfo->opposingFaceOrdinal_ << ")" << std::endl;  
      NaluEnv::self().naluOutput() << "  encapsulated by Gid: (";  
      for ( int ni = 0; ni < opposing_face_num_nodes; ++ni ) {
        stk::mesh::Entity node = opposing_face_node_rels[ni];
        NaluEnv::self().naluOutput() << bulk_data.identifier(node) << " ";
      }
      NaluEnv::self().naluOutput() << ")" << std::endl;
      NaluEnv::self().naluOutput() << "  INTERNAL CHECK.... does current Gp and opposing found Gp match coordiantes? What error?" << std::endl;
      NaluEnv::self().naluOutput() << "  current and opposing Gp coordinates:        " << std::endl;
      for ( int i = 0; i < nDim; ++i )
        NaluEnv::self().naluOutput() << "      " << i << " " << checkCurrentFaceGaussPointCoords[i] << " " << checkOpposingFaceGaussPointCoords[i] << std::endl;
      NaluEnv::self().naluOutput() << "  current and opposing Gp isoPar coordinates: " << std::endl;
      for ( int i = 0; i < nDim-1; ++i )
        NaluEnv::self().naluOutput() << "      " << i << " " << currentIsoParCoords[i] << " " << opposingIsoParCoords[i] << std::endl;
      NaluEnv::self().naluOutput() << std::endl;
      NaluEnv::self().naluOutput() << " in the end, the Error Distance Norm is: " << distanceNorm << std::endl;
      NaluEnv::self().naluOutput() << "-------------------------------------------------------------------" << std::endl;
    }
  }
}

} // namespace nalu
} // namespace sierra
