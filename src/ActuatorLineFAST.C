/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ActuatorLineFAST.h>
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
// ActuatorLineFASTInfo - holds all points in the tower specification
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineFASTInfo::ActuatorLineFASTInfo() 
  : processorId_(0),
    numPoints_(1),
    turbineName_("machine_one")
{
  // nothing to do
}


//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineFASTInfo::~ActuatorLineFASTInfo()
{
  // nothing to do
}


//==========================================================================
// Class Definition
//==========================================================================
// ActuatorLineFASTPointInfo - holds individual points information
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineFASTPointInfo::ActuatorLineFASTPointInfo( 
  size_t localId, 
  Point centroidCoords, 
  double searchRadius,
  Coordinates epsilon,
  double *velocity,
  ActuatorNodeType nType)
  : localId_(localId),
    centroidCoords_(centroidCoords),
    searchRadius_(searchRadius),
    epsilon_(epsilon),
    bestX_(1.0e16),
    bestElem_(stk::mesh::Entity()),
    nodeType_(nType)
{
  // initialize point velocity and displacement
  velocity_[0] = velocity[0];
  velocity_[1] = velocity[1];
  velocity_[2] = velocity[2];
}


//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineFASTPointInfo::~ActuatorLineFASTPointInfo()
{
  // nothing to do
}


//==========================================================================
// Class Definition
//==========================================================================
// ActuatorLineFAST - assemble source term for subgrin turbine; WIP
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineFAST::ActuatorLineFAST(
  Realm &realm,
  const YAML::Node &node)
  : Actuator(realm, node),
    realm_(realm),
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

     2) There can be many specifications with the number of points and omega processed.

     3) in the end, we fill the map of ActuatorLineFASTPointInfo objects and iterate this guy
        to assemble source terms

     4) at present, fake source terms on simple Gaussian weighting

    actuator:
      search_method: stk_octree
      search_target_part: block_1

      specifications:

        - name: machine_zero
          procNo: 0
          epsilon: [ 2.0, 0.0, 0.0 ]
          turbine_pos: [ 0.0, 0.0, 0.0 ]
          restart_filename: "blah"
          FAST_input_filename: "Test01.fst"
          turb_id:  1
          turbine_name: machine_zero
  */
}


//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineFAST::~ActuatorLineFAST()
{

  FAST.end(); // Call destructors in FAST_cInterface

  // delete data probes specifications vector
  for ( size_t k = 0; k < actuatorLineInfo_.size(); ++k )
    delete actuatorLineInfo_[k];
}


//--------------------------------------------------------------------------
//-------- compute_elem_force_given_weight ----------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::compute_elem_force_given_weight(
  const int &nDim,
  const double &g,
  const double *pointForce,
  double *elemForce)
{
  // Multiply the point force by the weight at this element location.
  for ( int j = 0; j < nDim; ++j )
    elemForce[j] = pointForce[j]*g;
}


//--------------------------------------------------------------------------
//-------- isotropic_Gaussian_projection -----------------------------------
//--------------------------------------------------------------------------
double
ActuatorLineFAST::isotropic_Gaussian_projection(
  const int &nDim,
  const double &dis,
  const Coordinates &epsilon)
{
  // Compute the force projection weight at this location using an
  // isotropic Gaussian.
  double g;
  const double pi = acos(-1.0);
  if ( nDim == 2 )
    g = (1.0 / (pow(epsilon.x_,2.0) * pi)) * exp(-pow((dis/epsilon.x_),2.0));
  else
    g = (1.0 / (pow(epsilon.x_,3.0) * pow(pi,1.5))) * exp(-pow((dis/epsilon.x_),2.0));

//std::cout << "g = " << g << std::endl;

  return g;
}


//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::load(
  const YAML::Node & y_node)
{
  // check for any data probes
  const YAML::Node y_actuatorLine = y_node["actuator"];
  if (y_actuatorLine) {
    NaluEnv::self().naluOutputP0() << "ActuatorLineFAST::load" << std::endl;

    // search specifications
    std::string searchMethodName = "na";
    get_if_present(y_actuatorLine, "search_method", searchMethodName, searchMethodName);
    
    // determine search method for this pair
    if ( searchMethodName == "boost_rtree" )
      searchMethod_ = stk::search::BOOST_RTREE;
    else if ( searchMethodName == "stk_octree" )
      searchMethod_ = stk::search::OCTREE;
    else if ( searchMethodName == "stk_kdtree" )
      searchMethod_ = stk::search::KDTREE;
    else
      NaluEnv::self().naluOutputP0() << "ActuatorLineFAST::search method not declared; will use BOOST_RTREE" << std::endl;

    // extract the set of from target names; each spec is homogeneous in this respect
    const YAML::Node searchTargets = y_actuatorLine["search_target_part"];
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

    try {
      FAST.readInputFile(y_actuatorLine);
    }
    catch( const std::runtime_error & ex) {
    }

    // save off number of towers
    const int nTurbinesGlob = FAST.get_nTurbinesGlob() ;

    // each specification can have multiple machines
    for (size_t iTurb = 0; iTurb < nTurbinesGlob; ++iTurb) {
      const YAML::Node cur_turbine = y_actuatorLine["Turbine"+std::to_string(iTurb)];
      ActuatorLineFASTInfo *actuatorLineInfo = new ActuatorLineFASTInfo();
      actuatorLineInfo_.push_back(actuatorLineInfo);
      
      // name
      const YAML::Node theName = cur_turbine["turbine_name"];
      if ( theName )
	actuatorLineInfo->turbineName_ = theName.as<std::string>() ;
      else
	throw std::runtime_error("ActuatorLineFAST: no name provided");
      
      // processor id - Get from FAST
      actuatorLineInfo->processorId_ = FAST.get_procNo(iTurb) ;
     
      actuatorLineMotion_ = true;
      
      // Force projection function properties
      const YAML::Node epsilon = cur_turbine["epsilon"];
        if ( epsilon )
          actuatorLineInfo->epsilon_ = epsilon.as<Coordinates>() ;
        else
          throw std::runtime_error("ActuatorLineFAST: lacking epsilon vector");
    }
  }
}


//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::setup()
{
  // objective: declare the part, register coordinates; must be before populate_mesh()

  double tStart;
  double tEnd;
  double dt ; 
  FAST.setRestart( realm_.restarted_simulation() ) ;
  tStart = realm_.get_current_time() ;
  FAST.setTstart( tStart ) ;
  dt = realm_.get_time_step_from_file();
  FAST.setDt ( dt ) ;
  if ( realm_.get_is_terminate_based_on_time() ) {
    tEnd = realm_.get_total_sim_time() ;
  } 
  else {
    const int ntEnd = realm_.get_max_time_step_count();
    const int ntStart = realm_.get_time_step_count() ;
    tEnd = (ntEnd - ntStart)*dt + tStart ;
  }
  FAST.setTend( tEnd );
  
  if ( ! FAST.isDryRun() ) {
    FAST.init() ;
  }

}


//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::initialize()
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  stk::mesh::MetaData & metaData = realm_.meta_data();
 
  const int nDim = metaData.spatial_dimension();

  // initialize need to ghost and elems to ghost
  needToGhostCount_ = 0;
  elemsToGhost_.clear();
  // initialize need to ghost and elems to ghost
  needToGhostForceCount_ = 0;
  elemsToGhostForce_.clear();

  // clear actuatorLinePointInfoMap_
  std::map<size_t, ActuatorLineFASTPointInfo *>::iterator iterPoint;
  for( iterPoint=actuatorLinePointInfoMap_.begin(); iterPoint!=actuatorLinePointInfoMap_.end(); ++iterPoint )
    delete (*iterPoint).second;
  actuatorLinePointInfoMap_.clear();
  for( iterPoint=actuatorLineForcePointInfoMap_.begin(); iterPoint!=actuatorLineForcePointInfoMap_.end(); ++iterPoint )
    delete (*iterPoint).second;
  actuatorLineForcePointInfoMap_.clear();
  
  bulkData.modification_begin();
  
  if ( actuatorLineGhosting_ == NULL) {
    // create new ghosting
    std::string theGhostName = "nalu_actuator_line_ghosting";
    actuatorLineGhosting_ = &bulkData.create_ghosting( theGhostName );
  }
  else {
    bulkData.destroy_ghosting(*actuatorLineGhosting_);
  }
  
  if ( actuatorLineForceGhosting_ == NULL) {
    // create new ghosting
    std::string theGhostName = "nalu_actuator_line_force_ghosting";
    actuatorLineForceGhosting_ = &bulkData.create_ghosting( theGhostName );
  }
  else {
    bulkData.destroy_ghosting(*actuatorLineForceGhosting_);
  }

  bulkData.modification_end();
  
  // clear some of the search info
  boundingSphereVec_.clear();
  boundingSphereForceVec_.clear();
  boundingElementBoxVec_.clear();
  searchKeyPair_.clear();
  searchKeyPairForce_.clear();

  // set all of the candidate elements in the search target names
  populate_candidate_elements();
  
  // create the ActuatorLineFASTPointInfo
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
ActuatorLineFAST::execute()
{
  // meta/bulk data and nDim
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const int nDim = metaData.spatial_dimension();

  // extract fields
  VectorFieldType *coordinates 
    = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  VectorFieldType *velocity = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  VectorFieldType *actuator_source 
    = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "actuator_source");
  ScalarFieldType *actuator_source_lhs
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "actuator_source_lhs");
  ScalarFieldType *g
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "g");
  ScalarFieldType *density
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density"); 
  // deal with proper viscosity
  const std::string viscName = realm_.is_turbulent() ? "effective_viscosity" : "viscosity";
  ScalarFieldType *viscosity
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, viscName); 

  // fixed size scratch
  std::vector<double> ws_pointGasVelocity(nDim);
  std::vector<double> ws_elemCentroid(nDim);
  std::vector<double> ws_pointForce(nDim);
  std::vector<double> ws_elemForce(nDim);
  double ws_pointGasDensity;
  double ws_pointGasViscosity;
  double ws_pointForceLHS;
  
  // zero out source term; do this manually since there are custom ghosted entities
  stk::mesh::Selector s_nodes = stk::mesh::selectField(*actuator_source);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * actSrc = stk::mesh::field_data(*actuator_source, b);
    double * actSrcLhs = stk::mesh::field_data(*actuator_source_lhs, b);
    double * gF = stk::mesh::field_data(*g, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      actSrcLhs[k] = 0.0;
      gF[k] = 0.0;
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

  // loop over map and get velocity at points
  std::map<size_t, ActuatorLineFASTPointInfo *>::iterator iterPoint;
  int np=0;
  for (iterPoint  = actuatorLinePointInfoMap_.begin();
       iterPoint != actuatorLinePointInfoMap_.end();
       ++iterPoint) {

    // actuator line info object of interest
    ActuatorLineFASTPointInfo * infoObject = (*iterPoint).second;

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

    if (FAST.isDebug() ) {
      NaluEnv::self().naluOutput() << "Node " << np << " Velocity = " << ws_pointGasVelocity[0] << " " << ws_pointGasVelocity[1] << " " << ws_pointGasVelocity[2] << " " << std::endl ;
    }

    FAST.setVelocity(ws_pointGasVelocity, np);
    np = np + 1;

  }    

  if ( ! FAST.isDryRun() ) {

    if ( FAST.isTimeZero() ) {
      FAST.solution0();
    }

    //Step FAST
    FAST.step();
  }
 
  // do we have mesh motion?
  if ( actuatorLineMotion_ )
    initialize();

  // loop over map and assemble source terms
  np = 0;
  for (iterPoint  = actuatorLineForcePointInfoMap_.begin();
       iterPoint != actuatorLineForcePointInfoMap_.end();
       ++iterPoint) {

    // actuator line info object of interest
    ActuatorLineFASTPointInfo * infoObject = (*iterPoint).second;

    //==========================================================================
    // extract the best element; compute drag given this velocity, property, etc
    // this point drag value will be used by all other elements below
    //==========================================================================
    stk::mesh::Entity bestElem = infoObject->bestElem_;
    int nodesPerElement = bulkData.num_nodes(bestElem);
    // compute volume
    double elemVolume = compute_volume(nDim, bestElem, bulkData);

    FAST.getForce(ws_pointForce, np);
    if (FAST.isDebug() ) {
      NaluEnv::self().naluOutput() << "Node " << np << " Type " << infoObject->nodeType_ << " Force = " << ws_pointForce[0] << " " << ws_pointForce[1] << " " << ws_pointForce[2] << " " << std::endl ;
    }
  //ws_pointForce[0] = 0.0; //Setting to zero for now
  //ws_pointForce[1] = 0.0;
  //ws_pointForce[2] = 0.0;
  //ws_pointForceLHS = 0.0; //Not clear what this should be - FIGURE IT OUT
    // assemble nodal quantity; radius should be zero, so we can apply fill point drag
  //assemble_source_to_nodes(nDim, bestElem, bulkData, elemVolume, &ws_pointForce[0], ws_pointForceLHS, 
  //                         *actuator_source, *actuator_source_lhs, 1.0);

    // get the vector of elements
    std::vector<stk::mesh::Entity> elementVec = infoObject->elementVec_;

    // Set up the necessary variables to check that forces/projection function are integrating up correctly.
    double gSum = 0.0;
    double forceSum[nDim];
    forceSum[0] = 0.0;
    forceSum[1] = 0.0;
    if (nDim>2){
      forceSum[2] = 0.0;
    }

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

      // compute distance
      const double distance = compute_distance(nDim, &ws_elemCentroid[0], &(infoObject->centroidCoords_[0]));

      double gA = 0.0;
      
      switch (infoObject->nodeType_) {
      case HUB:
        // project the force to this element centroid with projection function
        gA = isotropic_Gaussian_projection(nDim, distance, infoObject->epsilon_);
        compute_elem_force_given_weight(nDim, gA, &ws_pointForce[0], &ws_elemForce[0]);
        // assemble nodal quantity; no LHS contribution here...
        assemble_source_to_nodes(nDim, elem, bulkData, elemVolume, &ws_elemForce[0], ws_pointForceLHS,
                                 0.0, *actuator_source, *actuator_source_lhs, *g, 0.0);
        forceSum[0] += ws_elemForce[0]*elemVolume;
        forceSum[1] += ws_elemForce[1]*elemVolume;
        if (nDim > 2){
          forceSum[2] += ws_elemForce[2]*elemVolume;
        }
        gSum += gA*elemVolume;
        break;

      case BLADE:
        // project the force to this element centroid with projection function
        gA = isotropic_Gaussian_projection(nDim, distance, infoObject->epsilon_);
        compute_elem_force_given_weight(nDim, gA, &ws_pointForce[0], &ws_elemForce[0]);
        // assemble nodal quantity; no LHS contribution here...
        assemble_source_to_nodes(nDim, elem, bulkData, elemVolume, &ws_elemForce[0], ws_pointForceLHS,
                                 gA, *actuator_source, *actuator_source_lhs, *g, 0.0);
        forceSum[0] += ws_elemForce[0]*elemVolume;
        forceSum[1] += ws_elemForce[1]*elemVolume;
        if (nDim > 2){
          forceSum[2] += ws_elemForce[2]*elemVolume;
        }
        gSum += gA*elemVolume;
        break;

      case TOWER:
        // project the force to this element centroid with projection function
        gA = isotropic_Gaussian_projection(nDim, distance, infoObject->epsilon_);
        compute_elem_force_given_weight(nDim, gA, &ws_pointForce[0], &ws_elemForce[0]);
        // assemble nodal quantity; no LHS contribution here...
        assemble_source_to_nodes(nDim, elem, bulkData, elemVolume, &ws_elemForce[0], ws_pointForceLHS,
                                 gA, *actuator_source, *actuator_source_lhs, *g, 0.0);
        forceSum[0] += ws_elemForce[0]*elemVolume;
        forceSum[1] += ws_elemForce[1]*elemVolume;
        if (nDim > 2){
          forceSum[2] += ws_elemForce[2]*elemVolume;
        }
        gSum += gA*elemVolume;
        break;	
      }

    }

    NaluEnv::self().naluOutput() << "Actuator Point " << np << ", " << "# elems " << elementVec.size() << ", " << " Type " << infoObject->nodeType_ << std::endl;
    NaluEnv::self().naluOutput() << "  -Body force = " << forceSum[0] << " " << forceSum[1] << " " << forceSum[2] << " " << std::endl;
    NaluEnv::self().naluOutput() << "  -Act. force = " << ws_pointForce[0] << " " << ws_pointForce[1] << " " << ws_pointForce[2] << std::endl;
    NaluEnv::self().naluOutput() << "  -Ratio = " << forceSum[0]/ws_pointForce[0] << " " << forceSum[1]/ws_pointForce[1] << " " << forceSum[2]/ws_pointForce[2] << std::endl;
    NaluEnv::self().naluOutput() << "  -gSum = " << gSum << std::endl;

    np=np+1;
  }

  // parallel assemble (contributions from ghosted and locally owned)
  const std::vector<const stk::mesh::FieldBase*> sumFieldVec(1, actuator_source);
  stk::mesh::parallel_sum_including_ghosts(bulkData, sumFieldVec);

  const std::vector<const stk::mesh::FieldBase*> sumFieldG(1, g);
  stk::mesh::parallel_sum_including_ghosts(bulkData, sumFieldG);
}


//--------------------------------------------------------------------------
//-------- populate_candidate_elements -------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::populate_candidate_elements() 
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
      throw std::runtime_error("ActuatorLineFAST: Part is null" + searchTargetNames_[k]);     
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
ActuatorLineFAST::determine_elems_to_ghost()
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  stk::search::coarse_search(boundingSphereVec_, boundingElementBoxVec_, searchMethod_, 
                             NaluEnv::self().parallel_comm(), searchKeyPair_);
  stk::search::coarse_search(boundingSphereForceVec_, boundingElementBoxVec_, searchMethod_, 
                             NaluEnv::self().parallel_comm(), searchKeyPairForce_);

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

  for( ii=searchKeyPairForce_.begin(); ii!=searchKeyPairForce_.end(); ++ii ) {

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
      needToGhostForceCount_++;

      // deal with elements to push back to be ghosted
      stk::mesh::EntityProc theElemPair(theElemMeshObj, pt_proc);
      elemsToGhostForce_.push_back(theElemPair);
    }
  }

}

//--------------------------------------------------------------------------
//-------- create_actuator_line_point_info_map -----------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::create_actuator_line_point_info_map() {

  const double currentTime = realm_.get_current_time();

  stk::mesh::MetaData & metaData = realm_.meta_data(); 
  const int nDim = metaData.spatial_dimension();

  for ( size_t k = 0; k < actuatorLineInfo_.size(); ++k ) {
    
    const ActuatorLineFASTInfo *actuatorLineInfo = actuatorLineInfo_[k];
    
    int processorId = actuatorLineInfo->processorId_;
    if ( processorId == NaluEnv::self().parallel_rank() ) {
      
      // define a point that will hold the centroid
      Point centroidCoords;
    
      // local turbine id
      const int iTurbLoc = FAST.get_localTurbNo(k) ; // 'k' is the global turbine id here
      // scratch array for coordinates and dummy array for velocity
      double velocity[3] = {};
      double currentCoords[3] = {};

      // loop over all points for this turbine
      const int nBlades = FAST.get_numBlades(iTurbLoc); // Number of blades per turbine - only Turbine 0 for now
      const int numVelPtsBlade = FAST.get_numVelPtsBlade(iTurbLoc) ; // Number of nodes per blade - only Turbine 0 for now
      const int numVelPtsTwr = FAST.get_numVelPtsTwr(iTurbLoc) ; // Number of tower elements - only Turbine 0 for now
      const int numVelPts = FAST.get_numVelPts(iTurbLoc); // Total number of elements - only Turbine 0 for now
      const int numForcePtsBlade = FAST.get_numForcePtsBlade(iTurbLoc) ; // Number of nodes per blade - only Turbine 0 for now
      const int numForcePtsTwr = FAST.get_numForcePtsTwr(iTurbLoc) ; // Number of tower elements - only Turbine 0 for now
      const int numForcePts = FAST.get_numForcePts(iTurbLoc); // Total number of elements - only Turbine 0 for now
      //      actuatorLineInfo->numPoints_ = numVelPts ; // Can't change a const pointer
      
      if (! FAST.isDryRun() ) {
	int np = 0;
	for(int iNode = 0; iNode < numVelPts; iNode++) {
	  // extract current localPointId; increment for next one up...
	  size_t localPointId = localPointId_++;
	  stk::search::IdentProc<uint64_t,int> theIdent(localPointId, NaluEnv::self().parallel_rank());
	  
	  // set model coordinates from FAST
	  // move the coordinates; set the velocity... may be better on the lineInfo object
	  FAST.getVelNodeCoordinates(currentCoords, iNode);
	  if (FAST.isDebug() ) {
	    NaluEnv::self().naluOutput() << "Vel Node " << np << " Position = " << currentCoords[0] << " " << currentCoords[1] << " " << currentCoords[2] << " " << std::endl ;
	  }
	  velocity[0] = 0.0;
	  velocity[1] = 0.0;
	  velocity[2] = 0.0;

          double searchRadius = actuatorLineInfo->epsilon_.x_ * sqrt(log(1.0/0.001));
	  
	  for ( int j = 0; j < nDim; ++j )
	    centroidCoords[j] = currentCoords[j];
	  
	  // create the bounding point sphere and push back
	  boundingSphere theSphere( Sphere(centroidCoords, searchRadius), theIdent);
	  boundingSphereVec_.push_back(theSphere);
	  
	  // create the point info and push back to map
	  ActuatorLineFASTPointInfo *actuatorLinePointInfo 
	    = new ActuatorLineFASTPointInfo(localPointId, centroidCoords, 
					searchRadius, actuatorLineInfo->epsilon_,
                                        velocity, FAST.getVelNodeType(k, iNode));
	  actuatorLinePointInfoMap_[localPointId] = actuatorLinePointInfo;
	  
	  np=np+1;
	}

	np = 0;
	for(int iNode = 0; iNode < numForcePts; iNode++) {
	  // extract current localPointId; increment for next one up...
	  size_t localPointId = localPointId_++;
	  stk::search::IdentProc<uint64_t,int> theIdent(localPointId, NaluEnv::self().parallel_rank());
	  
	  // set model coordinates from FAST
	  // move the coordinates; set the velocity... may be better on the lineInfo object
	  FAST.getForceNodeCoordinates(currentCoords, iNode);
	  if (FAST.isDebug() ) {
	    NaluEnv::self().naluOutput() << "Force Node " << np << " Position = " << currentCoords[0] << " " << currentCoords[1] << " " << currentCoords[2] << " " << std::endl ;
	  }
	  velocity[0] = 0.0;
	  velocity[1] = 0.0;
	  velocity[2] = 0.0;

          double searchRadius = actuatorLineInfo->epsilon_.x_ * sqrt(log(1.0/0.001));
	  
	  for ( int j = 0; j < nDim; ++j )
	    centroidCoords[j] = currentCoords[j];
	  
	  // create the bounding point sphere and push back
	  boundingSphere theSphere( Sphere(centroidCoords, searchRadius), theIdent);
	  boundingSphereForceVec_.push_back(theSphere);
	  
	  // create the point info and push back to map
	  ActuatorLineFASTPointInfo *actuatorLinePointInfo 
	    = new ActuatorLineFASTPointInfo(localPointId, centroidCoords, 
					searchRadius, actuatorLineInfo->epsilon_,
                                        velocity, FAST.getForceNodeType(k, iNode));
	  actuatorLineForcePointInfoMap_[localPointId] = actuatorLinePointInfo;
	  
	  np=np+1;
	}

      }
      else {
	NaluEnv::self().naluOutput() << "Proc " << NaluEnv::self().parallel_rank() << " loc iTurb " << iTurbLoc << " glob iTurb " << k << std::endl ;
      }
      
    }

  }
}


//--------------------------------------------------------------------------
//-------- manage_ghosting -------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::manage_ghosting() 
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  
  // check for ghosting need
  uint64_t g_needToGhostCount = 0;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &needToGhostCount_, &g_needToGhostCount, 1);
  if (g_needToGhostCount > 0) {
    NaluEnv::self().naluOutputP0() << "ActuatorLineFAST alg will ghost a number of entities for velocity actuator nodes: "
                                   << g_needToGhostCount  << std::endl;
    bulkData.modification_begin();
    bulkData.change_ghosting( *actuatorLineGhosting_, elemsToGhost_);
    bulkData.modification_end();
  }
  else {
    NaluEnv::self().naluOutputP0() << "ActuatorLineFAST alg will NOT ghost entities for velcity actuator nodes: " << std::endl;
  }

  uint64_t g_needToGhostForceCount = 0;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &needToGhostForceCount_, &g_needToGhostForceCount, 1);
  if (g_needToGhostForceCount > 0) {
    NaluEnv::self().naluOutputP0() << "ActuatorLineFAST alg will ghost a number of entities for force actuator nodes: "
                                   << g_needToGhostForceCount  << std::endl;
    bulkData.modification_begin();
    bulkData.change_ghosting( *actuatorLineForceGhosting_, elemsToGhostForce_);
    bulkData.modification_end();
  }
  else {
    NaluEnv::self().naluOutputP0() << "ActuatorLineFAST alg will NOT ghost entities for force actuator nodes: " << std::endl;
  }

}
 
 
//--------------------------------------------------------------------------
//-------- complete_search -------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::complete_search()
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
      std::map<size_t, ActuatorLineFASTPointInfo *>::iterator iterPoint;
      iterPoint=actuatorLinePointInfoMap_.find(thePt);
      if ( iterPoint == actuatorLinePointInfoMap_.end() )
        throw std::runtime_error("no valid entry for actuatorLinePointInfoMap_");
      
      // extract the point object and push back the element to either the best 
      // candidate or the standard vector of elements
      ActuatorLineFASTPointInfo *actuatorLinePointInfo = iterPoint->second;

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
        actuatorLinePointInfo->bestElem_ = elem;
      }
        actuatorLinePointInfo->elementVec_.push_back(elem);
    }
    else {
      // not this proc's issue
    }
  }

  for( ii=searchKeyPairForce_.begin(); ii!=searchKeyPairForce_.end(); ++ii ) {

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
      std::map<size_t, ActuatorLineFASTPointInfo *>::iterator iterPoint;
      iterPoint=actuatorLineForcePointInfoMap_.find(thePt);
      if ( iterPoint == actuatorLineForcePointInfoMap_.end() )
        throw std::runtime_error("no valid entry for actuatorLineForcePointInfoMap_");
      
      // extract the point object and push back the element to either the best 
      // candidate or the standard vector of elements
      ActuatorLineFASTPointInfo *actuatorLineForcePointInfo = iterPoint->second;

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
                                                        &(actuatorLineForcePointInfo->centroidCoords_[0]),
                                                        &(isoParCoords[0]));
 

      // save off best element and its isoparametric coordinates for this point
      if ( nearestDistance < actuatorLineForcePointInfo->bestX_ ) {
        actuatorLineForcePointInfo->bestX_ = nearestDistance;
        actuatorLineForcePointInfo->isoParCoords_ = isoParCoords;
        actuatorLineForcePointInfo->bestElem_ = elem;
      }
      actuatorLineForcePointInfo->elementVec_.push_back(elem);
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
ActuatorLineFAST::resize_std_vector( 
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
ActuatorLineFAST::gather_field(
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
ActuatorLineFAST::gather_field_for_interp(
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
ActuatorLineFAST::compute_volume(
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
ActuatorLineFAST::interpolate_field(
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
ActuatorLineFAST::compute_elem_centroid(
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
//-------- compute_distance ------------------------------------------------
//--------------------------------------------------------------------------
double
ActuatorLineFAST::compute_distance(
  const int &nDim,
  const double *elemCentroid,
  const double *pointCentroid) 
{ 
  double distance = 0.0;
  for ( int j = 0; j < nDim; ++j )
    distance += std::pow(elemCentroid[j] - pointCentroid[j], 2);
  distance = std::sqrt(distance);
  return distance;
}


//--------------------------------------------------------------------------
//-------- assemble_source_to_nodes ----------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineFAST::assemble_source_to_nodes(
  const int &nDim,
  stk::mesh::Entity elem,
  const stk::mesh::BulkData & bulkData,
  const double &elemVolume,
  const double *drag,
  const double &dragLHS,
  const double &gLocal,
  stk::mesh::FieldBase &actuator_source,
  stk::mesh::FieldBase &actuator_source_lhs,
  stk::mesh::FieldBase &g,
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
    double * sourceTerm = (double*)stk::mesh::field_data(actuator_source, node );
    double * sourceTermLHS = (double*)stk::mesh::field_data(actuator_source_lhs, node );
    double * gGlobal = (double*)stk::mesh::field_data(g, node);


    // nodal weight based on volume weight
    const double nodalWeight = ws_scv_volume_[ip]/elemVolume;
    *sourceTermLHS += nodalWeight*dragLHS*lhsFac;
    *gGlobal += gLocal;
    for ( int j=0; j < nDim; ++j ) {
      sourceTerm[j] += nodalWeight*drag[j];
    }
  }
} 
 

} // namespace nalu
} // namespace Sierra
