/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ActuatorLine.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>
#include <NaluEnv.h>
#include <Realm.h>
#include <Simulation.h>

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
// ActuatorLineInfo - holds all points in the tower specification
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineInfo::ActuatorLineInfo() 
  : processorId_(0),
    numPoints_(1),
    turbineName_("machine_one"),
    radius_(0),
    omega_(0.0),
    twoSigSq_(0.0)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineInfo::~ActuatorLineInfo()
{
  // nothing to do
}

//==========================================================================
// Class Definition
//==========================================================================
// ActuatorLinePointInfo - holds individual points information
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLinePointInfo::ActuatorLinePointInfo( 
  size_t localId, 
  Point centroidCoords, 
  double radius, 
  double omega,
  double twoSigSq) 
  : localId_(localId),
    centroidCoords_(centroidCoords),
    radius_(radius),
    omega_(omega),
    twoSigSq_(twoSigSq),
    bestX_(1.0e16),
    bestElem_(stk::mesh::Entity())
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLinePointInfo::~ActuatorLinePointInfo()
{
  // nothing to do
}

//==========================================================================
// Class Definition
//==========================================================================
// ActuatorLine - assemble source term for subgrin turbine; WIP
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLine::ActuatorLine(
  Realm &realm,
  const YAML::Node &node)
  : realm_(realm),
    searchMethod_(stk::search::BOOST_RTREE),
    actuatorLineGhosting_(NULL),
    needToGhostCount_(0),
    localPointId_(0),
    actuatorLineMotion_(false),
    pi_(acos(-1.0))
{
  // load the data
  load(node);

  /*
    current WIP prototype
    Design concepts:
     1) First and foremost, elements are ghosted to the owning point rank. This 
        probably should be changed since the number of elements might be larger
        than the number of points. Therefore, ghosting points to elements is probably
        easier. This will remove the parallel sum contributions from ghosted elements.
        time will tell..

     2) There can be many specifications with the number of points and omaga processed.

     3) in the end, we fill the map of ActuatorLinePointInfo objects and iterate this guy
        to assemble source terms

     4) at present, fake source terms on simple gaussian weighting

    actuator_line:
      search_method: stk_octree
      search_target_part: block_1

      specifications:

        - name: machine_one 
          radius: 2.0
          omega: 1.0
          gaussian_decay_radius: 1.5
          gaussian_decay_target: 0.01
          coodinates: [0.0, 0.0, 0.0]
  */
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLine::~ActuatorLine()
{
  // delete data probes specifications vector
  for ( size_t k = 0; k < actuatorLineInfo_.size(); ++k )
    delete actuatorLineInfo_[k];
}

//--------------------------------------------------------------------------
//-------- compute_point_drag ----------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::compute_point_drag( 
  const int &nDim,
  const double &pointRadius, 
  const double *pointGasVelocity,
  const double &pointGasViscosity,
  const double &pointGasDensity,
  double *pointDrag,
  double &pointDragLHS)
{
  // HACK... assume point velocity is zero
  double pointVelocity = 0.0;
  double vRelMag = 0.0;
  for ( int j = 0; j < nDim; ++j )
    vRelMag += (pointVelocity - pointGasVelocity[j] )*(pointVelocity - pointGasVelocity[j]);
  vRelMag = std::sqrt(vRelMag);

  // Reynolds number and friction factors
  double ReP = 2.0*pointGasDensity*pointRadius*vRelMag/pointGasViscosity;
  double CubeRtReP = (ReP < 1000.) ? std::cbrt(ReP) : 0.0;
  double fD = (ReP < 1000.0) ? (1.0 + CubeRtReP*CubeRtReP/6.0) : (0.424/24.0 * ReP);
  double coef = 6.0*pi_*pointGasViscosity*pointRadius;

  // this is from the fluids perspective, not the psuego particle
  pointDragLHS = 2.0*coef*fD;
  for ( int j = 0; j < nDim; ++j )
    pointDrag[j] = coef*fD*(pointVelocity - pointGasVelocity[j]);
}

//--------------------------------------------------------------------------
//-------- compute_elem_drag_given_radius ----------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::compute_elem_drag_given_radius( 
  const int &nDim,
  const double &radius, 
  const double &twoSigSq,
  const double *pointDrag,
  double *elemDrag)
{
  // gaussian weight based on radius
  const double gaussWeight = std::exp(-radius*radius/twoSigSq);
  for ( int j = 0; j < nDim; ++j )
    elemDrag[j] = pointDrag[j]*gaussWeight;
}

//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::load(
  const YAML::Node & y_node)
{
  // check for any data probes
  const YAML::Node *y_actuatorLine = y_node.FindValue("actuator_line");
  if (y_actuatorLine) {
    NaluEnv::self().naluOutputP0() << "ActuatorLine::load" << std::endl;

    // search specifications
    std::string searchMethodName = "na";
    get_if_present(*y_actuatorLine, "search_method", searchMethodName, searchMethodName);
    
    // determine search method for this pair
    if ( searchMethodName == "boost_rtree" )
      searchMethod_ = stk::search::BOOST_RTREE;
    else if ( searchMethodName == "stk_octree" )
      searchMethod_ = stk::search::OCTREE;
    else
      NaluEnv::self().naluOutputP0() << "ActuatorLine::search method not declared; will use BOOST_RTREE" << std::endl;

    // extract the set of from target names; each spec is homogeneous in this respect
    const YAML::Node &searchTargets = (*y_actuatorLine)["search_target_part"];
    if (searchTargets.Type() == YAML::NodeType::Scalar) {
      searchTargetNames_.resize(1);
      searchTargets >> searchTargetNames_[0];
    }
    else {
      searchTargetNames_.resize(searchTargets.size());
      for (size_t i=0; i < searchTargets.size(); ++i) {
        searchTargets[i] >> searchTargetNames_[i];
      }
    }
    
    const YAML::Node *y_specs = expect_sequence(*y_actuatorLine, "specifications", false);
    if (y_specs) {

      // save off number of towers
      const int numTowers = y_specs->size();

      // deal with processors... Distribute each tower over subsequent procs
      const int numProcs = NaluEnv::self().parallel_size();
      const int divProcTower = std::max(numProcs/numTowers, numProcs);

      // each specification can have multiple machines
      for (size_t ispec = 0; ispec < y_specs->size(); ++ispec) {
        const YAML::Node &y_spec = (*y_specs)[ispec];
        
        ActuatorLineInfo *actuatorLineInfo = new ActuatorLineInfo();
        actuatorLineInfo_.push_back(actuatorLineInfo);
        
        // name
        const YAML::Node *theName = y_spec.FindValue("turbine_name");
        if ( theName )
          *theName >> actuatorLineInfo->turbineName_;
        else
          throw std::runtime_error("ActuatorLine: no name provided");
        
        // processor id; distribute los equally over the number of processors
        actuatorLineInfo->processorId_ = divProcTower > 0 ? ispec % divProcTower : 0;

        // number of points
        get_if_present(y_spec, "number_of_points", actuatorLineInfo->numPoints_, actuatorLineInfo->numPoints_);
        if ( actuatorLineInfo->numPoints_ > 1 )
          throw std::runtime_error("ActuatorLine: number of points must be unity");

        // radius and omega
        get_if_present(y_spec, "radius", actuatorLineInfo->radius_, actuatorLineInfo->radius_);
        get_if_present(y_spec, "omega", actuatorLineInfo->omega_, actuatorLineInfo->omega_);
        if ( actuatorLineInfo->omega_ != 0.0 )
          throw std::runtime_error("ActuatorLine: not ready for omega not equal to zero");
        
        // finally, the gaussian props
        double gaussDecayTarget = 0.01;
        double gaussDecayRadius = 1.5;
        get_if_present(y_spec, "gaussian_decay_radius", gaussDecayRadius, gaussDecayRadius);
        get_if_present(y_spec, "gaussian_decay_target", gaussDecayTarget, gaussDecayTarget);
        actuatorLineInfo->twoSigSq_ = -gaussDecayRadius*gaussDecayRadius/std::log(gaussDecayTarget);
        
        // coordinates of this point
        const YAML::Node *coord = y_spec.FindValue("coordinates");
        if ( coord )
          *coord >> actuatorLineInfo->coordinates_;
        else
          throw std::runtime_error("ActuatorLine: lacking coordinates");
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::setup()
{
  // objective: declare the part, register coordinates; must be before populate_mesh()
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::initialize()
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  stk::mesh::MetaData & metaData = realm_.meta_data();
 
  const int nDim = metaData.spatial_dimension();

  // initialize need to ghost and elems to ghost
  needToGhostCount_ = 0;
  elemsToGhost_.clear();

  // clear actuatorLinePointInfoMap_
  std::map<size_t, ActuatorLinePointInfo *>::iterator iterPoint;
  for( iterPoint=actuatorLinePointInfoMap_.begin(); iterPoint!=actuatorLinePointInfoMap_.end(); ++iterPoint )
    delete (*iterPoint).second;
  actuatorLinePointInfoMap_.clear();
  
  bulkData.modification_begin();
  
  if ( actuatorLineGhosting_ == NULL) {
    // create new ghosting
    std::string theGhostName = "nalu_actuator_line_ghosting";
    actuatorLineGhosting_ = &bulkData.create_ghosting( theGhostName );
  }
  else {
    bulkData.destroy_ghosting(*actuatorLineGhosting_);
  }
  
  bulkData.modification_end();
  
  // clear some of the search info
  boundingSphereVec_.clear();
  boundingElementBoxVec_.clear();
  searchKeyPair_.clear();

  // set all of the candidate elements in the search target names
  populate_candidate_elements();
  
  // create the ActuatorLinePointInfo
  create_actuator_line_point_info_map();

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
ActuatorLine::execute()
{
  // do we have mesh motion?
  if ( actuatorLineMotion_ )
    initialize();

  // meta/bulk data and nDim
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const int nDim = metaData.spatial_dimension();

  // extract fields
  VectorFieldType *coordinates 
    = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  VectorFieldType *velocity = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  VectorFieldType *actuator_line_source 
    = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "actuator_line_source");
  ScalarFieldType *actuator_line_source_lhs
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "actuator_line_source_lhs");
  ScalarFieldType *density
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density"); 
  // deal with proper viscosity
  const std::string viscName = realm_.is_turbulent() ? "effective_viscosity" : "viscosity";
  ScalarFieldType *viscosity
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, viscName); 

  // fixed size scratch
  std::vector<double> ws_pointGasVelocity(nDim);
  std::vector<double> ws_elemCentroid(nDim);
  std::vector<double> ws_pointDrag(nDim);
  std::vector<double> ws_elemDrag(nDim);
  double ws_pointGasDensity;
  double ws_pointGasViscosity;
  double ws_pointDragLHS;
  
  // zero out source term; do this manually since there are custom ghosted entities
  stk::mesh::Selector s_nodes = stk::mesh::selectField(*actuator_line_source);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * actSrc = stk::mesh::field_data(*actuator_line_source, b);
    double * actSrcLhs = stk::mesh::field_data(*actuator_line_source_lhs, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      actSrcLhs[k] = 0.0;
      const int offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        actSrc[offSet+j] = 0.0;
      }
    }
  }

  // parallel communicate data to the ghosted elements; again can communicate points to element ranks
  if ( NULL != actuatorLineGhosting_ ) {
    std::vector< const stk::mesh::FieldBase *> ghostFieldVec;
    // fields that are needed
    ghostFieldVec.push_back(coordinates);
    ghostFieldVec.push_back(velocity);
    ghostFieldVec.push_back(viscosity);
    stk::mesh::communicate_field_data(*actuatorLineGhosting_, ghostFieldVec);
  }

  // loop over map and assemble source terms
  std::map<size_t, ActuatorLinePointInfo *>::iterator iterPoint;
  for (iterPoint  = actuatorLinePointInfoMap_.begin();
       iterPoint != actuatorLinePointInfoMap_.end();
       ++iterPoint) {

    // actuator line info object of interest
    ActuatorLinePointInfo * infoObject = (*iterPoint).second;

    //==========================================================================
    // extract the best element; compute drag given this velocity, property, etc
    // this point drag value will be used by all other elements below
    //==========================================================================
    stk::mesh::Entity bestElem = infoObject->bestElem_;
    int nodesPerElement = bulkData.num_nodes(bestElem);
    
    // resize some work vectors
    resize_std_vector(nDim, ws_coordinates_, bestElem, bulkData);
    resize_std_vector(nDim, ws_velocity_, bestElem, bulkData);
    resize_std_vector(1, ws_viscosity_, bestElem, bulkData);
    resize_std_vector(1, ws_density_, bestElem, bulkData);

    // gather nodal data to element nodes; both vector and scalar; coords are used in determinant calc
    gather_field(nDim, &ws_coordinates_[0], *coordinates, bulkData.begin_nodes(bestElem), 
                 nodesPerElement);
    gather_field_for_interp(nDim, &ws_velocity_[0], *velocity, bulkData.begin_nodes(bestElem), 
                            nodesPerElement);
    gather_field_for_interp(1, &ws_viscosity_[0], *viscosity, bulkData.begin_nodes(bestElem), 
                            nodesPerElement);
    gather_field_for_interp(1, &ws_density_[0], *density, bulkData.begin_nodes(bestElem), 
                            nodesPerElement);

    // compute volume
    double elemVolume = compute_volume(nDim, bestElem, bulkData);

    // interpolate velocity
    interpolate_field(nDim, bestElem, bulkData, &(infoObject->isoParCoords_[0]), 
                      &ws_velocity_[0], &ws_pointGasVelocity[0]);
    
    // interpolate viscosity
    interpolate_field(1, bestElem, bulkData, &(infoObject->isoParCoords_[0]), 
                      &ws_viscosity_[0], &ws_pointGasViscosity);

    // interpolate density
    interpolate_field(1, bestElem, bulkData, &(infoObject->isoParCoords_[0]), 
                      &ws_density_[0], &ws_pointGasDensity);
    
    // point drag calculation
    compute_point_drag(nDim, infoObject->radius_, &ws_pointGasVelocity[0], ws_pointGasViscosity, 
                       ws_pointGasDensity, &ws_pointDrag[0], ws_pointDragLHS);
        
    // assemble nodal quantity; radius should be zero, so we can apply fill point drag
    assemble_source_to_nodes(nDim, bestElem, bulkData, elemVolume, &ws_pointDrag[0], ws_pointDragLHS, 
                             *actuator_line_source, *actuator_line_source_lhs, 1.0);

    // get the vector of elements
    std::vector<stk::mesh::Entity> elementVec = infoObject->elementVec_;

    // iterate them and apply source term; gather coords
    for ( size_t k = 0; k < elementVec.size(); ++k ) {

      stk::mesh::Entity elem = elementVec[k];

      nodesPerElement = bulkData.num_nodes(elem);
    
      // resize some work vectors
      resize_std_vector(nDim, ws_coordinates_, elem, bulkData);

      // gather coordinates
      gather_field(nDim, &ws_coordinates_[0], *coordinates, bulkData.begin_nodes(elem), 
                   nodesPerElement);

      // compute volume
      double elemVolume = compute_volume(nDim, elem, bulkData);

      // determine element centroid
      compute_elem_centroid(nDim, &ws_elemCentroid[0], nodesPerElement);

      // compute radius
      const double radius = compute_radius(nDim, &ws_elemCentroid[0], &(infoObject->centroidCoords_[0]));
    
      // get drag at this element centroid with proper Gaussian weighting
      compute_elem_drag_given_radius(nDim, radius, infoObject->twoSigSq_, &ws_pointDrag[0], &ws_elemDrag[0]);
    
      // assemble nodal quantity; no LHS contribution here...
      assemble_source_to_nodes(nDim, elem, bulkData, elemVolume, &ws_elemDrag[0], ws_pointDragLHS, 
                               *actuator_line_source, *actuator_line_source_lhs, 0.0);
    } 
  }

  // parallel assemble (contributions from ghosted and locally owned)
  const std::vector<const stk::mesh::FieldBase*> sumFieldVec(1, actuator_line_source);
  stk::mesh::parallel_sum_including_ghosts(bulkData, sumFieldVec);
}

//--------------------------------------------------------------------------
//-------- populate_candidate_elements -------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::populate_candidate_elements() 
{
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  const int nDim = metaData.spatial_dimension();

  // fields
  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // point data structures
  Point minCorner, maxCorner;

  // extract part
  stk::mesh::PartVector searchParts;
  for ( size_t k = 0; k < searchTargetNames_.size(); ++k ) {
    stk::mesh::Part *thePart = metaData.get_part(searchTargetNames_[k]);
    if ( NULL != thePart )
      searchParts.push_back(thePart);
    else
      throw std::runtime_error("ActuatorLine: Part is null" + searchTargetNames_[k]);     
  }

  // selector and bucket loop
  stk::mesh::Selector s_locally_owned = metaData.locally_owned_part()
    &stk::mesh::selectUnion(searchParts);
  
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
ActuatorLine::determine_elems_to_ghost()
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  stk::search::coarse_search(boundingSphereVec_, boundingElementBoxVec_, searchMethod_, 
                             NaluEnv::self().parallel_comm(), searchKeyPair_);

  // lowest effort is to ghost elements to the owning rank of the point; can just as easily do the opposite
  std::vector<std::pair<boundingSphere::second_type, boundingElementBox::second_type> >::const_iterator ii;
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
//-------- create_actuator_line_point_info_map -----------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::create_actuator_line_point_info_map() {
  stk::mesh::MetaData & metaData = realm_.meta_data(); 
  const int nDim = metaData.spatial_dimension();

  for ( size_t k = 0; k < actuatorLineInfo_.size(); ++k ) {
    
    const ActuatorLineInfo *actuatorLineInfo = actuatorLineInfo_[k];
    
    int processorId = actuatorLineInfo->processorId_;
    if ( processorId == NaluEnv::self().parallel_rank() ) {
      
      // define a point that will hold the centroid
      Point centroidCoords;
      
      // loop over all points
      for ( int j = 0; j < actuatorLineInfo->numPoints_; ++j ) {
        // extract current localPointId; increment for next one up...
        size_t localPointId = localPointId_++;
        stk::search::IdentProc<uint64_t,int> theIdent(localPointId, NaluEnv::self().parallel_rank());
        
        // extract model coordinates
        centroidCoords[0] = actuatorLineInfo->coordinates_.x_;
        centroidCoords[1] = actuatorLineInfo->coordinates_.y_;
        if ( nDim > 2 )
          centroidCoords[2] = actuatorLineInfo->coordinates_.z_;
        
        // move the coordinates; set the velocity... may be better on the lineInfo object
        set_current_coordinates(actuatorLineInfo->omega_);
        set_current_velocity(actuatorLineInfo->omega_);

        // create the bounding point sphere and push back
        boundingSphere theSphere( Sphere(centroidCoords, actuatorLineInfo->radius_), theIdent);
        boundingSphereVec_.push_back(theSphere);
        
        // create the point info and push back to map
        ActuatorLinePointInfo *actuatorLinePointInfo 
          = new ActuatorLinePointInfo(localPointId, centroidCoords, 
                                      actuatorLineInfo->radius_, actuatorLineInfo->omega_, 
                                      actuatorLineInfo->twoSigSq_);
        actuatorLinePointInfoMap_[localPointId] = actuatorLinePointInfo;
      }
    }
  }  
}

//--------------------------------------------------------------------------
//-------- set_current_coordinates -----------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::set_current_coordinates() 
{
  // to do
}

//--------------------------------------------------------------------------
//-------- set_current_coordinates -----------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::set_current_coordinates( const double &/*omega*/) 
{
  // to do
}

//--------------------------------------------------------------------------
//-------- set_current_velocity --------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::set_current_velocity( const double &/*omega*/) 
{
  // to do
}

//--------------------------------------------------------------------------
//-------- manage_ghosting -------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::manage_ghosting() 
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  
  // check for ghosting need
  uint64_t g_needToGhostCount = 0;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &needToGhostCount_, &g_needToGhostCount, 1);
  if (g_needToGhostCount > 0) {
    NaluEnv::self().naluOutputP0() << "ActuatorLine alg will ghost a number of entities: "
                                   << g_needToGhostCount  << std::endl;
    bulkData.modification_begin();
    bulkData.change_ghosting( *actuatorLineGhosting_, elemsToGhost_);
    bulkData.modification_end();
  }
  else {
    NaluEnv::self().naluOutputP0() << "ActuatorLine alg will NOT ghost entities: " << std::endl;
  }
}
  
//--------------------------------------------------------------------------
//-------- complete_search -------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::complete_search()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const int nDim = metaData.spatial_dimension();

  // extract fields
  VectorFieldType *coordinates 
    = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
 
  // now proceed with the standard search
  std::vector<std::pair<boundingSphere::second_type, boundingElementBox::second_type> >::const_iterator ii;
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
      std::map<size_t, ActuatorLinePointInfo *>::iterator iterPoint;
      iterPoint=actuatorLinePointInfoMap_.find(thePt);
      if ( iterPoint == actuatorLinePointInfoMap_.end() )
        throw std::runtime_error("no valid entry for actuatorLinePointInfoMap_");
      
      // extract the point object and push back the element to either the best 
      // candidate or the standard vector of elements
      ActuatorLinePointInfo *actuatorLinePointInfo = iterPoint->second;

      // extract topo and master element for this topo
      const stk::mesh::Bucket &theBucket = bulkData.bucket(elem);
      const stk::topology &elemTopo = theBucket.topology();
      MasterElement *meSCS = realm_.get_surface_master_element(elemTopo);
      const int nodesPerElement = meSCS->nodesPerElement_;

      // gather elemental coords
      std::vector<double> elementCoords(nDim*nodesPerElement);
      gather_field(nDim, &elementCoords[0], *coordinates, bulkData.begin_nodes(elem), 
                   nodesPerElement);

      
      // find isoparametric points
      std::vector<double> isoParCoords(nDim);
      const double nearestDistance = meSCS->isInElement(&elementCoords[0],
                                                        &(actuatorLinePointInfo->centroidCoords_[0]),
                                                        &(isoParCoords[0]));
      
      // save off best element and its isoparametric coordinates for this point
      if ( nearestDistance < actuatorLinePointInfo->bestX_ ) { 
        actuatorLinePointInfo->bestX_ = nearestDistance;    
        actuatorLinePointInfo->isoParCoords_ = isoParCoords;
        if ( stk::mesh::Entity() == actuatorLinePointInfo->bestElem_ ) {
          actuatorLinePointInfo->bestElem_ = elem;
        }
        else { 
          // swap the current best element
          actuatorLinePointInfo->elementVec_.push_back(actuatorLinePointInfo->bestElem_);
          actuatorLinePointInfo->bestElem_ = elem;
        }
      }   
      else {
        // regular element
        actuatorLinePointInfo->elementVec_.push_back(elem);
      }
    }
    else {
      // not this proc's issue
    }
  }
}

//--------------------------------------------------------------------------
//-------- resize_std_vector -----------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::resize_std_vector( 
  const int &sizeOfField,
  std::vector<double> &theVector,   
  stk::mesh::Entity elem, 
  const stk::mesh::BulkData & bulkData)
{
  const stk::topology &elemTopo = bulkData.bucket(elem).topology();
  MasterElement *meSCS = realm_.get_surface_master_element(elemTopo);
  const int nodesPerElement = meSCS->nodesPerElement_;
  theVector.resize(nodesPerElement*sizeOfField);
}

//--------------------------------------------------------------------------
//-------- gather_field ----------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::gather_field(
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
//-------- gather_field_for_interp -----------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::gather_field_for_interp(
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
      const int offSet = j*nodesPerElement + ni; 
      fieldToFill[offSet] = theField[j];
    }   
  }   
}

//--------------------------------------------------------------------------
//-------- compute_volume --------------------------------------------------
//--------------------------------------------------------------------------
double
ActuatorLine::compute_volume(
  const int &nDim,
  stk::mesh::Entity elem, 
  const stk::mesh::BulkData & bulkData) 
{
  // extract master element from the bucket in which the element resides
  const stk::topology &elemTopo = bulkData.bucket(elem).topology();
  MasterElement *meSCV = realm_.get_volume_master_element(elemTopo);
  int nodesPerElement = meSCV->nodesPerElement_;
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
//-------- interpolate_field -----------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::interpolate_field(
  const int &sizeOfField,
  stk::mesh::Entity elem, 
  const stk::mesh::BulkData & bulkData,
  double *isoParCoords,
  const double *fieldAtNodes,
  double *pointField) 
{
  // extract master element from the bucket in which the element resides
  const stk::topology &elemTopo = bulkData.bucket(elem).topology();
  MasterElement *meSCS = realm_.get_surface_master_element(elemTopo);
  
  // interpolate velocity to this best point
  meSCS->interpolatePoint(
    sizeOfField,
    isoParCoords,
    fieldAtNodes,
    pointField); 
}

//--------------------------------------------------------------------------
//-------- compute_elem_centroid -------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLine::compute_elem_centroid(
  const int &nDim,
  double *elemCentroid,
  const int & nodesPerElement) 
{
  // zero
  for ( int j = 0; j < nDim; ++j )
    elemCentroid[j] = 0.0;
  
  // assemble
  for ( int ni = 0; ni < nodesPerElement; ++ni ) {
    for ( int j=0; j < nDim; ++j ) {
      elemCentroid[j] += ws_coordinates_[ni*nDim+j]/nodesPerElement;
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_radius ------------------------------------------------
//--------------------------------------------------------------------------
double
ActuatorLine::compute_radius(
  const int &nDim,
  const double *elemCentroid,
  const double *pointCentroid) 
{ 
  double radius = 0.0;
  for ( int j = 0; j < nDim; ++j )
    radius += std::pow(elemCentroid[j] - pointCentroid[j], 2);
  radius = std::sqrt(radius);
  return radius;
}

//--------------------------------------------------------------------------
//-------- assemble_source_to_nodes ----------------------------------------
//-------------------------------------------------------------------------
void
ActuatorLine::assemble_source_to_nodes(
  const int &nDim,
  stk::mesh::Entity elem, 
  const stk::mesh::BulkData & bulkData,
  const double &elemVolume,
  const double *drag,
  const double &dragLHS,
  stk::mesh::FieldBase &actuator_line_source,
  stk::mesh::FieldBase &actuator_line_source_lhs,
  const double &lhsFac) 
{
  // extract master element from the bucket in which the element resides
  const stk::topology &elemTopo = bulkData.bucket(elem).topology();
  MasterElement *meSCV = realm_.get_volume_master_element(elemTopo);
  int nodesPerElement = meSCV->nodesPerElement_;
  const int numScvIp = meSCV->numIntPoints_;

  // extract elem_node_relations
  stk::mesh::Entity const* elem_node_rels = bulkData.begin_nodes(elem);

  // assemble to nodes
  const int *ipNodeMap = meSCV->ipNodeMap();
  for ( int ip = 0; ip < numScvIp; ++ip ) {
      
    // nearest node to ip
    const int nearestNode = ipNodeMap[ip];
    
    // extract node and pointer to source term
    stk::mesh::Entity node = elem_node_rels[nearestNode];
    double * sourceTerm = (double*)stk::mesh::field_data(actuator_line_source, node );
    double * sourceTermLHS = (double*)stk::mesh::field_data(actuator_line_source_lhs, node );
    
    // nodal weight based on volume weight
    const double nodalWeight = ws_scv_volume_[ip]/elemVolume;
    *sourceTermLHS += nodalWeight*dragLHS*lhsFac;
    for ( int j=0; j < nDim; ++j ) {
      sourceTerm[j] += nodalWeight*drag[j];
    }
  } 
}
  
} // namespace nalu
} // namespace Sierra
