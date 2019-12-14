/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <ComputeWallFrictionVelocityProjectedAlgorithm.h>
#include <Algorithm.h>
#include <PointInfo.h>

#include <FieldTypeDef.h>
#include <Realm.h>
#include <master_element/MasterElement.h>
#include <NaluEnv.h>

#include <utils/StkHelpers.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/ExodusTranslator.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

// basic c++
#include <cmath>

namespace sierra{
namespace nalu{

// compare operator - move this....
struct compareId {
  bool operator () (const std::pair<uint64IdentProc, uint64IdentProc> &p, const uint64_t i) {
    return (p.first.id() < i);
  }
  bool operator () (const uint64_t i, const std::pair<uint64IdentProc, uint64IdentProc> &p) {
    return (i < p.first.id());
  }
};

struct lessThan
{
  bool operator() (const std::pair<uint64IdentProc, uint64IdentProc> &p, const std::pair<uint64IdentProc, uint64IdentProc> &q) {
    return (p.first.id() < q.first.id());
  }
};

//==========================================================================
// Class Definition
//==========================================================================
// ComputeWallFrictionVelocityProjectedAlgorithm - utau at wall bc
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ComputeWallFrictionVelocityProjectedAlgorithm::ComputeWallFrictionVelocityProjectedAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  const double projectedDistance,
  const bool useShifted,
  std::vector<std::vector<PointInfo *> > &pointInfoVec,
  stk::mesh::Ghosting *wallFunctionGhosting)
  : Algorithm(realm, part),
    useShifted_(useShifted),
    pointInfoVec_(pointInfoVec),
    wallFunctionGhosting_(wallFunctionGhosting),
    bulkData_(&realm.bulk_data()),
    metaData_(&realm.meta_data()),
    nDim_(realm.meta_data().spatial_dimension()),
    yplusCrit_(11.63),
    elog_(9.8),
    kappa_(realm.get_turb_model_constant(TM_kappa)),
    maxIteration_(20),
    tolerance_(1.0e-6),
    firstInitialization_(true),
    provideOutput_(false),
    searchMethod_(stk::search::KDTREE),
    expandBoxPercentage_(0.05),
    needToGhostCount_(0)
{
  // save off fields
  velocity_ = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  bcVelocity_ = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "wall_velocity_bc");
  coordinates_ = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  density_ = metaData_->get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  viscosity_ = metaData_->get_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity");
  exposedAreaVec_ = metaData_->get_field<GenericFieldType>(metaData_->side_rank(), "exposed_area_vector");
  wallFrictionVelocityBip_ = metaData_->get_field<GenericFieldType>(metaData_->side_rank(), "wall_friction_velocity_bip");
  wallNormalDistanceBip_ = metaData_->get_field<GenericFieldType>(metaData_->side_rank(), "wall_normal_distance_bip");
  assembledWallArea_ = metaData_->get_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_wall_area_wf");
  assembledWallNormalDistance_ = metaData_->get_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_wall_normal_distance");
  
  // set data
  set_data(projectedDistance);

  // what do we need ghosted for this alg to work?
  ghostFieldVec_.push_back(&(velocity_->field_of_state(stk::mesh::StateNP1)));
}
  
//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ComputeWallFrictionVelocityProjectedAlgorithm::~ComputeWallFrictionVelocityProjectedAlgorithm()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeWallFrictionVelocityProjectedAlgorithm::execute()
{
  // fixed size
  std::vector<double> uProjected(nDim_);
  std::vector<double> uBcBip(nDim_);
  std::vector<double> unitNormal(nDim_);
  std::vector<double> cProjected(nDim_);
  
  // pointers to fixed values
  double *p_uProjected = &uProjected[0];
  double *p_uBcBip = &uBcBip[0];
  double *p_unitNormal= &unitNormal[0];
  
  // isopar coordinates for the owning element
  std::vector<double> isoParCoords(nDim_);
  
  // nodal fields to gather
  std::vector<double> ws_bcVelocity;
  std::vector<double> ws_density;
  std::vector<double> ws_viscosity;
  
  // master element
  std::vector<double> ws_face_shape_function;
  
  // deal with state
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);
  
  // define vector of parent topos; should always be UNITY in size
  std::vector<stk::topology> parentTopo;

  // simple for now...
  initialize();
  
  // parallel communicate ghosted entities
  if ( nullptr != wallFunctionGhosting_ )
    stk::mesh::communicate_field_data(*(wallFunctionGhosting_), ghostFieldVec_);
  
  // iterate over parts to match construction (requires global counter over locally owned faces)
  size_t pointInfoVecCounter = 0;  
  for ( size_t pv = 0; pv < partVec_.size(); ++pv ) {

    // extract projected distance
    const double pDistance = projectedDistanceVec_[pv];

    // define selector (per part)
    stk::mesh::Selector s_locally_owned 
      = metaData_->locally_owned_part() &stk::mesh::Selector(*partVec_[pv]);
    
    stk::mesh::BucketVector const& face_buckets =
      realm_.get_buckets( metaData_->side_rank(), s_locally_owned );
    
    for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
          ib != face_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      
      // face master element
      MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());
      const int nodesPerFace = b.topology().num_nodes();
      const int numScsBip = meFC->numIntPoints_;
      
      // mapping from ip to nodes for this ordinal; face perspective (use with face_node_relations)
      const int *faceIpNodeMap = meFC->ipNodeMap();

      // algorithm related; element
      ws_bcVelocity.resize(nodesPerFace*nDim_);
      ws_density.resize(nodesPerFace);
      ws_viscosity.resize(nodesPerFace);
      ws_face_shape_function.resize(numScsBip*nodesPerFace);
      
      // pointers
      double *p_bcVelocity = &ws_bcVelocity[0];
      double *p_density = &ws_density[0];
      double *p_viscosity = &ws_viscosity[0];
      double *p_face_shape_function = &ws_face_shape_function[0];
      
      // shape functions
      if ( useShifted_ )
        meFC->shifted_shape_fcn(&p_face_shape_function[0]);
      else
        meFC->shape_fcn(&p_face_shape_function[0]);

      const stk::mesh::Bucket::size_type length   = b.size();
      
      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        
        // get face
        stk::mesh::Entity face = b[k];
        
        //======================================
        // gather nodal data off of face
        //======================================
        stk::mesh::Entity const * face_node_rels = bulkData_->begin_nodes(face);
        int num_face_nodes = bulkData_->num_nodes(face);
        // sanity check on num nodes
        ThrowAssert( num_face_nodes == nodesPerFace );
        for ( int ni = 0; ni < num_face_nodes; ++ni ) {
          stk::mesh::Entity node = face_node_rels[ni];
          
          // gather scalars
          p_density[ni]    = *stk::mesh::field_data(densityNp1, node);
          p_viscosity[ni] = *stk::mesh::field_data(*viscosity_, node);
          
          // gather vectors
          double * uBc = stk::mesh::field_data(*bcVelocity_, node);
          const int niNdim = ni*nDim_;
          for ( int j=0; j < nDim_; ++j ) {
            p_bcVelocity[niNdim+j] = uBc[j];
          }
        }
        
        // pointer to face data
        const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, face);
        double *wallNormalDistanceBip = stk::mesh::field_data(*wallNormalDistanceBip_, face);
        double *wallFrictionVelocityBip = stk::mesh::field_data(*wallFrictionVelocityBip_, face);
        
        // extract the vector of PointInfo for this face 
        std::vector<PointInfo *> &faceInfoVec = pointInfoVec_[pointInfoVecCounter++];

        // loop over ips
        for ( int ip = 0; ip < numScsBip; ++ip ) {

          // extract nearest node to this ip
          const int localFaceNode = faceIpNodeMap[ip]; 
          stk::mesh::Entity nearestNode = face_node_rels[localFaceNode];
          
          // extract point info for this ip - must matches the construction of the pInfo vector
          PointInfo *pInfo = faceInfoVec[ip];
          stk::mesh::Entity owningElement = pInfo->owningElement_;

          // get master element type for this contactInfo
          MasterElement *meSCS  = pInfo->meSCS_;
          const int nodesPerElement = meSCS->nodesPerElement_;
          std::vector <double > elemNodalVelocity(nodesPerElement*nDim_);
          std::vector <double > elemNodalCoords(nodesPerElement*nDim_);
          std::vector <double > shpfc(nodesPerElement);

          // gather element data
          stk::mesh::Entity const* elem_node_rels = bulkData_->begin_nodes(owningElement);
          const int num_elem_nodes = bulkData_->num_nodes(owningElement);
          for ( int ni = 0; ni < num_elem_nodes; ++ni ) {
            stk::mesh::Entity node = elem_node_rels[ni];
            // gather velocity (conforms to interpolatePoint)
            const double *uNp1 = stk::mesh::field_data(velocityNp1, node );
            const double *coords = stk::mesh::field_data(*coordinates_, node );
            for ( int j = 0; j < nDim_; ++j ) {
              elemNodalVelocity[j*nodesPerElement+ni] = uNp1[j];
              elemNodalCoords[j*nodesPerElement+ni] = coords[j];
            }
          }
          
          // interpolate to elemental point location
          meSCS->interpolatePoint(
            nDim_,
            &(pInfo->isoParCoords_[0]),
            &elemNodalVelocity[0],
            &uProjected[0]);

          // sanity check for coords
          meSCS->interpolatePoint(
            nDim_,
            &(pInfo->isoParCoords_[0]),
            &elemNodalCoords[0],
            &cProjected[0]);

          if ( provideOutput_ ) {
            for (int j = 0; j < nDim_; ++j ) 
              NaluEnv::self().naluOutput() << "Coords sanity check: " << cProjected[j] << " " << pInfo->pointCoordinates_[j] << std::endl;
          }
          
          // zero out vector quantities; squeeze in aMag
          double aMag = 0.0;
          for ( int j = 0; j < nDim_; ++j ) {
            p_uBcBip[j] = 0.0;
            const double axj = areaVec[ip*nDim_+j];
            aMag += axj*axj;
          }
          aMag = std::sqrt(aMag);
          
          // interpolate to bip
          double rhoBip = 0.0;
          double muBip = 0.0;
          const int ipNpf = ip*nodesPerFace;
          for ( int ic = 0; ic < nodesPerFace; ++ic ) {
            const double r = p_face_shape_function[ipNpf+ic];
            rhoBip += r*p_density[ic];
            muBip += r*p_viscosity[ic];
            const int icNdim = ic*nDim_;
            for ( int j = 0; j < nDim_; ++j ) {
              p_uBcBip[j] += r*p_bcVelocity[icNdim+j];
            }
          }

          // form unit normal and determine yp (approximated by 1/4 distance along edge)
          for ( int j = 0; j < nDim_; ++j ) {
            p_unitNormal[j] = areaVec[ip*nDim_+j]/aMag;
          }
          
          double ypBip = pDistance;
          wallNormalDistanceBip[ip] = ypBip;
          
          // assemble to nodal quantities
          double * assembledWallArea = stk::mesh::field_data(*assembledWallArea_, nearestNode );
          double * assembledWallNormalDistance = stk::mesh::field_data(*assembledWallNormalDistance_, nearestNode );
          
          *assembledWallArea += aMag;
          *assembledWallNormalDistance += aMag*ypBip;

          // determine tangential velocity
          double uTangential = 0.0;
          for ( int i = 0; i < nDim_; ++i ) {
            double uiTan = 0.0;
            double uiBcTan = 0.0;
            for ( int j = 0; j < nDim_; ++j ) {
              const double ninj = p_unitNormal[i]*p_unitNormal[j];
              if ( i==j ) {
                const double om_nini = 1.0 - ninj;
                uiTan += om_nini*p_uProjected[j];
                uiBcTan += om_nini*p_uBcBip[j];
              }
              else {
                uiTan -= ninj*p_uProjected[j];
                uiBcTan -= ninj*p_uBcBip[j];
              }
            }
            uTangential += (uiTan-uiBcTan)*(uiTan-uiBcTan);
          }
          uTangential = std::sqrt(uTangential);
          
          // provide an initial guess based on yplusCrit_ (more robust than a pure guess on utau)
          double utauGuess = yplusCrit_*muBip/rhoBip/ypBip;
          
          compute_utau(uTangential, ypBip, rhoBip, muBip, utauGuess);
          wallFrictionVelocityBip[ip] = utauGuess;
          
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- set_data --------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeWallFrictionVelocityProjectedAlgorithm::set_data( 
  double theDouble)
{
  projectedDistanceVec_.push_back(theDouble);
}

//--------------------------------------------------------------------------
//-------- compute_utau ----------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeWallFrictionVelocityProjectedAlgorithm::compute_utau(
    const double &up, const double &yp,
    const double &density, const double &viscosity,
    double &utau )
{
  bool converged = false;

  const double A = elog_*density*yp/viscosity;

  for ( int k = 0; k < maxIteration_; ++k ) {

    const double wrk = std::log(A*utau);

    // evaluate F'

    const double fPrime = -(1.0+wrk);

    // evaluate function
    const double f = kappa_*up- utau*wrk;

    // update variable
    const double df = f/fPrime;

    utau -= df;
    if ( std::abs(df) < tolerance_ ) {
      converged = true;
      break;
    }
  }

  // report trouble
  if (!converged ) {
    NaluEnv::self().naluOutputP0() << "Issue with utau; not converged " << std::endl;
    NaluEnv::self().naluOutputP0() << up << " " << yp << " " << utau << std::endl;
  }
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeWallFrictionVelocityProjectedAlgorithm::initialize()
{
  
  // only process if first time or mesh motion is active
  if ( !firstInitialization_  && !realm_.has_mesh_motion() )
    return;
  
  // clear some of the search info
  boundingPointVec_.clear();
  boundingBoxVec_.clear();
  searchKeyPair_.clear();

  // initialize all ghosting data structures
  initialize_ghosting();
  
  // construct if the size is zero; reset always
  if ( pointInfoVec_.size() == 0 )
    construct_bounding_points();
  reset_point_info();
  
  // construct the bounding boxes
  construct_bounding_boxes();

  // coarse search (fills elemsToGhost_)
  coarse_search();

  manage_ghosting();

  // complete search
  complete_search();
  
  // set flag for the next possible time we are through the initialization method
  firstInitialization_ = false;
}

//--------------------------------------------------------------------------
//-------- initialize_ghosting ---------------------------------------------
//--------------------------------------------------------------------------
void
ComputeWallFrictionVelocityProjectedAlgorithm::initialize_ghosting()
{
  // initialize need to ghost and elems to ghost
  needToGhostCount_ = 0;
  elemsToGhost_.clear();
  
  bulkData_->modification_begin();  
  if ( nullptr == wallFunctionGhosting_) {
    // create new ghosting
    std::string theGhostName = "nalu_wall_function_ghosting";
    wallFunctionGhosting_ = &(bulkData_->create_ghosting( theGhostName ));
  }
  else {
    bulkData_->destroy_ghosting(*wallFunctionGhosting_);
  }
  bulkData_->modification_end();
}

//--------------------------------------------------------------------------
//-------- construct_bounding_points --------------------------------------
//--------------------------------------------------------------------------
void
ComputeWallFrictionVelocityProjectedAlgorithm::construct_bounding_points()
{
  // hold the point location projected off the face integration points
  Point ipCoordinates;
  Point pointCoordinates;

  // field extraction
  VectorFieldType *coordinates
    = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  
  // nodal fields to gather
  std::vector<double> ws_coordinates;
  std::vector<double> ws_face_shape_function;

  // fixed size
  std::vector<double> ws_unitNormal(nDim_,0.0);
  
  // need to keep track of some sort of local id for each gauss point...
  uint64_t localPointId = 0;
  
  // iterate over parts to allow for projected distance to vary per surface and defines ordering everywhere else
  for ( size_t pv = 0; pv < partVec_.size(); ++pv ) {
    
    // extract projected distance
    const double pDistance = projectedDistanceVec_[pv];
    
    // define selector (per part)
    stk::mesh::Selector s_locally_owned 
      = metaData_->locally_owned_part() &stk::mesh::Selector(*partVec_[pv]);
    
    stk::mesh::BucketVector const& face_buckets =
      realm_.get_buckets( metaData_->side_rank(), s_locally_owned );
    
    for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
          ib != face_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      
      // face master element
      MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());
      const int nodesPerFace = b.topology().num_nodes();
      const int numScsBip = meFC->numIntPoints_;
      
      // algorithm related; element
      ws_coordinates.resize(nodesPerFace*nDim_);
      ws_face_shape_function.resize(numScsBip*nodesPerFace);
      
      // pointers
      double *p_coordinates = &ws_coordinates[0];
      double *p_face_shape_function = &ws_face_shape_function[0];
      
      // shape functions
      if ( useShifted_ )
        meFC->shifted_shape_fcn(&p_face_shape_function[0]);
      else
        meFC->shape_fcn(&p_face_shape_function[0]);
      
      const stk::mesh::Bucket::size_type length   = b.size();
      
      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        
        // get face
        stk::mesh::Entity face = b[k];
        
        // pointer to face data
        const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, face);
        
        //======================================
        // gather nodal data off of face
        //======================================
        stk::mesh::Entity const * face_node_rels = bulkData_->begin_nodes(face);
        int num_face_nodes = bulkData_->num_nodes(face);
        // sanity check on num nodes
        ThrowAssert( num_face_nodes == nodesPerFace );
        for ( int ni = 0; ni < num_face_nodes; ++ni ) {
          stk::mesh::Entity node = face_node_rels[ni]; 
          // gather vectors
          double * coords = stk::mesh::field_data(*coordinates, node);
          const int niNdim = ni*nDim_;
          for ( int j=0; j < nDim_; ++j ) {
            p_coordinates[niNdim+j] = coords[j];
          }
        }
                        
        // set size for vector of points on this face
        std::vector<PointInfo *> faceInfoVec(numScsBip);
        for ( int ip = 0; ip < numScsBip; ++ip ) { 
          
          // compute area magnitude
          double aMag = 0.0;
          for ( int j = 0; j < nDim_; ++j ) {
            const double axj = areaVec[ip*nDim_+j];
            aMag += axj*axj;
          }
          aMag = std::sqrt(aMag);
          
          // compute normal (outward facing)
          for ( int j = 0; j < nDim_; ++j ) {
            ws_unitNormal[j] = areaVec[ip*nDim_+j]/aMag;
            ipCoordinates[j] = 0.0;
          }
          
          // interpolate coodinates to gauss point
          const int ipNpf = ip*nodesPerFace;
          for ( int ic = 0; ic < nodesPerFace; ++ic ) {
            const double r = p_face_shape_function[ipNpf+ic];
            for ( int j = 0; j < nDim_; ++j ) {
              ipCoordinates[j] += r*p_coordinates[ic*nDim_+j];
            }
          }
          
          // project in space (unit normal is outward facing, hence the -)
          for ( int j = 0; j < nDim_; ++j ) {
            pointCoordinates[j] = ipCoordinates[j] - pDistance*ws_unitNormal[j];
          }
          
          // setup ident for this point; use local integration point id
          uint64IdentProc theIdent(localPointId, NaluEnv::self().parallel_rank());
          
          // create the bounding point and push back
          boundingPoint bPoint(Point(pointCoordinates), theIdent);
          boundingPointVec_.push_back(bPoint);
          
          PointInfo *pInfo = new PointInfo(bPoint, localPointId, ipCoordinates, pointCoordinates, nDim_);
          faceInfoVec[ip] = pInfo;
          localPointId++;
        }
        
        // push them all back
        pointInfoVec_.push_back(faceInfoVec);
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- construct_bounding_boxes ----------------------------------------
//--------------------------------------------------------------------------
void
ComputeWallFrictionVelocityProjectedAlgorithm::construct_bounding_boxes()
{
  // extract coordinates
  VectorFieldType *coordinates
    = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  
  // setup data structures for search
  Point minCorner, maxCorner;

  // deal with element part vector
  stk::mesh::PartVector elemBlockPartVec;
  const stk::mesh::PartVector allParts = metaData_->get_parts();
  for ( size_t k = 0; k < allParts.size(); ++k ) {
    if ( stk::mesh::is_element_block(*allParts[k]) ) {
      elemBlockPartVec.push_back(allParts[k]);
    }
  }
  
  // selector
  stk::mesh::Selector s_locally_owned_union
    = metaData_->locally_owned_part() &stk::mesh::selectUnion(elemBlockPartVec);
  stk::mesh::BucketVector const &elem_buckets 
    = bulkData_->get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );
  
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib;
    
    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      // get element
      stk::mesh::Entity element = b[k];
      
      // initialize max and min
      for (int j = 0; j < nDim_; ++j ) {
        minCorner[j] = +1.0e16;
        maxCorner[j] = -1.0e16;
      }
      
      // extract elem_node_relations
      stk::mesh::Entity const* elem_node_rels = bulkData_->begin_nodes(element);
      const int num_nodes = bulkData_->num_nodes(element);
      
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];
        
        // pointers to real data
        const double * coords = stk::mesh::field_data(*coordinates, node );
        
        // check max/min
        for ( int j = 0; j < nDim_; ++j ) {
          minCorner[j] = std::min(minCorner[j], coords[j]);
          maxCorner[j] = std::max(maxCorner[j], coords[j]);
        }
      }
      
      // setup ident
      uint64IdentProc theIdent(bulkData_->identifier(element), NaluEnv::self().parallel_rank());
      
      // expand the box 
      for ( int i = 0; i < nDim_; ++i ) {
        const double theMin = minCorner[i];
        const double theMax = maxCorner[i];
        const double increment = expandBoxPercentage_*(theMax - theMin);
        minCorner[i] -= increment;
        maxCorner[i] += increment;
      }

      // correct for 2d
      //if ( nDim_ == 2 ) {
      //  minCorner[2] = -1.0;
      //  maxCorner[2] = +1.0;
      // }
      
      // create the bounding box and push back
      boundingBox theBox(Box(minCorner,maxCorner), theIdent);
      boundingBoxVec_.push_back(theBox);
    }
  }
}

//--------------------------------------------------------------------------
//-------- reset_point_info -----------------------------------------------
//--------------------------------------------------------------------------
void
ComputeWallFrictionVelocityProjectedAlgorithm::reset_point_info()
{
  std::vector<std::vector<PointInfo*> >::iterator ii;
  for( ii=pointInfoVec_.begin(); ii!=pointInfoVec_.end(); ++ii ) {
    std::vector<PointInfo *> &theVec = (*ii);    
    for ( size_t k = 0; k < theVec.size(); ++k ) {
      PointInfo *pInfo = theVec[k];
      pInfo->bestX_ = pInfo->bestXRef_;
    }
  }
}

//--------------------------------------------------------------------------
//-------- coarse_search ---------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeWallFrictionVelocityProjectedAlgorithm::coarse_search()
{
  // first coarse search
  stk::search::coarse_search(boundingPointVec_, boundingBoxVec_,
                             searchMethod_, NaluEnv::self().parallel_comm(), searchKeyPair_);
  
  // sort the product of the search
  std::sort(searchKeyPair_.begin(), searchKeyPair_.end(), lessThan());
  
  // now determine elements to ghost
  std::vector<std::pair<boundingPoint::second_type, boundingBox::second_type> >::const_iterator ii;  
  for ( ii=searchKeyPair_.begin(); ii!=searchKeyPair_.end(); ++ii ) {
    const uint64_t theBox = ii->second.id();
    unsigned theRank = NaluEnv::self().parallel_rank();
    const unsigned pt_proc = ii->first.proc();
    const unsigned box_proc = ii->second.proc();
    if ( (box_proc == theRank) && (pt_proc != theRank) ) {
      
      // find the element
      stk::mesh::Entity theElemMeshObj = bulkData_->get_entity(stk::topology::ELEMENT_RANK, theBox);
      if ( !(bulkData_->is_valid(theElemMeshObj)) )
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
//-------- manage_ghosting -------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeWallFrictionVelocityProjectedAlgorithm::manage_ghosting()
{  
  // check for ghosting need
  size_t g_needToGhostCount = 0;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &needToGhostCount_, &g_needToGhostCount, 1);
  if (g_needToGhostCount > 0) {
    
    NaluEnv::self().naluOutputP0() << "Projected LOW alg will ghost a number of entities: "
                                   << g_needToGhostCount  << std::endl;
    
    bulkData_->modification_begin();
    bulkData_->change_ghosting( *wallFunctionGhosting_, elemsToGhost_);
    bulkData_->modification_end();
    // no linsys contribution and, hence, no populate_ghost_comm_procs()
  }
  else {
    NaluEnv::self().naluOutputP0() << "Projected LOW alg will NOT ghost entities: " << std::endl;
  }
  
  // ensure that the coordinates for the ghosted elements (required for the fine search) are up-to-date
  if (g_needToGhostCount > 0 ) {
    VectorFieldType *coordinates 
      = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
    std::vector<const stk::mesh::FieldBase*> fieldVec = {coordinates};
    stk::mesh::communicate_field_data(*wallFunctionGhosting_, fieldVec);
  }
}

//--------------------------------------------------------------------------
//-------- complete_search--------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeWallFrictionVelocityProjectedAlgorithm::complete_search()
{
  // fields
  VectorFieldType *coordinates = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // coordinates
  std::vector<double> isoParCoords(nDim_);
  std::vector<double> pointCoords(nDim_);

  // invert the process... Loop over InfoVec_ and query searchKeyPair_ for this information (avoids a map)
  std::vector<PointInfo *> problemInfoVec;
  std::vector<std::vector<PointInfo *> >::iterator ii;
  for( ii=pointInfoVec_.begin(); ii!=pointInfoVec_.end(); ++ii ) {
    std::vector<PointInfo *> &theVec = (*ii);
    for ( size_t k = 0; k < theVec.size(); ++k ) {
      
      PointInfo *pInfo = theVec[k];
      const uint64_t localPointId  = pInfo->localPointId_; 
      for ( int j = 0; j < nDim_; ++j )
        pointCoords[j] = pInfo->pointCoordinates_[j];
      
      std::pair <std::vector<std::pair<uint64IdentProc, uint64IdentProc> >::const_iterator, std::vector<std::pair<uint64IdentProc, uint64IdentProc> >::const_iterator > 
        p2 = std::equal_range(searchKeyPair_.begin(), searchKeyPair_.end(), localPointId, compareId());
      
      if ( p2.first == p2.second ) {
        problemInfoVec.push_back(pInfo);        
      }
      else {
        for (std::vector<std::pair<uint64IdentProc, uint64IdentProc> >::const_iterator jj = p2.first; jj != p2.second; ++jj ) {
          
          const uint64_t theBox = jj->second.id();
          const unsigned theRank = NaluEnv::self().parallel_rank();
          const unsigned pt_proc = jj->first.proc();
        
          // check if I own the point...
          if ( theRank == pt_proc ) {

            // proceed as required; all elements should have already been ghosted via the coarse search
            stk::mesh::Entity candidateElement = bulkData_->get_entity(stk::topology::ELEMENT_RANK, theBox);
            if ( !(bulkData_->is_valid(candidateElement)) )
              throw std::runtime_error("no valid entry for element");
            
            int elemIsGhosted = bulkData_->bucket(candidateElement).owned() ? 0 : 1;
                        
            // now load the elemental nodal coords
            stk::mesh::Entity const * elem_node_rels = bulkData_->begin_nodes(candidateElement);
            int num_nodes = bulkData_->num_nodes(candidateElement);            
            std::vector<double> elementCoords(nDim_*num_nodes);
            
            for ( int ni = 0; ni < num_nodes; ++ni ) {
              stk::mesh::Entity node = elem_node_rels[ni];
              // gather coordinates (conforms to isInElement)
              const double * coords =  stk::mesh::field_data(*coordinates, node);
              for ( int j = 0; j < nDim_; ++j ) {
                elementCoords[j*num_nodes+ni] = coords[j];
              }
            }
            
            // extract the topo from this element...
            const stk::topology elemTopo = bulkData_->bucket(candidateElement).topology();
            MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(elemTopo);

            const double nearestDistance = meSCS->isInElement(&elementCoords[0],
                                                              &(pointCoords[0]),
                                                              &(isoParCoords[0]));

            // check if this element is the best
            if ( nearestDistance < pInfo->bestX_ ) {
              pInfo->owningElement_ = candidateElement;
              pInfo->meSCS_ = meSCS;
              pInfo->isoParCoords_ = isoParCoords;
              pInfo->bestX_ = nearestDistance;
              pInfo->elemIsGhosted_ = elemIsGhosted;
            }
          }
          else {
            // not this proc's issue
          }
        }
      }
    }
  }
  
  if ( provideOutput_ ) {
    
    // provide output
    for( ii=pointInfoVec_.begin(); ii!=pointInfoVec_.end(); ++ii ) {
      std::vector<PointInfo *> &theVec = (*ii);
      for ( size_t k = 0; k < theVec.size(); ++k ) {        
        provide_output(theVec[k], false);        
      }
    }
    
    // sanity check on the elements provided in the bounding box
    for ( size_t k = 0; k < boundingBoxVec_.size(); ++k ) {   
      NaluEnv::self().naluOutput() << "element ids provided to search: " << boundingBoxVec_[k].second.id() << std::endl; 
      Box bB = boundingBoxVec_[k].first;
      NaluEnv::self().naluOutput() <<  ".... x: " << bB.get_x_min() << " " << bB.get_x_max() << std::endl;
      NaluEnv::self().naluOutput() <<  ".... y: " << bB.get_y_min() << " " << bB.get_y_max() << std::endl;
      NaluEnv::self().naluOutput() <<  ".... z: " << bB.get_z_min() << " " << bB.get_z_max() << std::endl;
    }
  }

  // check if there was a problem
  if ( problemInfoVec.size() > 0 ) {
    NaluEnv::self().naluOutputP0() << "there was BIG PROBLEM with problemInfoVec " << std::endl; 
    
    for ( size_t k = 0; k < problemInfoVec.size(); ++k ) {   
      provide_output(problemInfoVec[k], true);        
    }
    
    // sanity check on the elements provided in the bounding box
    for ( size_t k = 0; k < boundingBoxVec_.size(); ++k ) {   
      NaluEnv::self().naluOutput() << " BIG PROBLEM: element ids provided to search: " << boundingBoxVec_[k].second.id() << std::endl; 
      Box bB = boundingBoxVec_[k].first;
      NaluEnv::self().naluOutput() <<  ".... x: " << bB.get_x_min() << " " << bB.get_x_max() << std::endl;
      NaluEnv::self().naluOutput() <<  ".... y: " << bB.get_y_min() << " " << bB.get_y_max() << std::endl;
      NaluEnv::self().naluOutput() <<  ".... z: " << bB.get_z_min() << " " << bB.get_z_max() << std::endl;
    }
    throw std::runtime_error("ComputeWallFrictionVelocityProjectedAlgorithm::search_error()");
  }
}

//--------------------------------------------------------------------------
//-------- provide_output --------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeWallFrictionVelocityProjectedAlgorithm::provide_output( 
  const PointInfo *pInfo,
  const bool problemPoint)
{
  const uint64_t localId = pInfo->localPointId_;
  stk::mesh::Entity theElem = pInfo->owningElement_;
  
  NaluEnv::self().naluOutput() << "...Review for Point ip: " << localId << std::endl;
  if ( problemPoint )
    NaluEnv::self().naluOutput() << "   BAD POINT" << std::endl;
  
  NaluEnv::self().naluOutput() << "   owning element id:    " << bulkData_->identifier(theElem) << std::endl;
  
  NaluEnv::self().naluOutput() << "   face coordinates: ";
  for ( int j = 0; j < nDim_; ++j )
    NaluEnv::self().naluOutput() << " " << pInfo->ipCoordinates_[j] << " ";
  NaluEnv::self().naluOutput() << std::endl;
  
  NaluEnv::self().naluOutput() << "   proj coordinates: ";
  for ( int j = 0; j < nDim_; ++j )
    NaluEnv::self().naluOutput() << " " << pInfo->pointCoordinates_[j] << " ";
  NaluEnv::self().naluOutput() << std::endl;
  
  NaluEnv::self().naluOutput() << "   point coordinates: " << std::endl;
  Point bP = pInfo->bPoint_.first;
  NaluEnv::self().naluOutput() <<  ".... x: " << bP.get_x_min() << " " << bP.get_x_max() << std::endl;
  NaluEnv::self().naluOutput() <<  ".... y: " << bP.get_y_min() << " " << bP.get_y_max() << std::endl;
  NaluEnv::self().naluOutput() <<  ".... z: " << bP.get_z_min() << " " << bP.get_z_max() << std::endl;   
}

} // namespace nalu
} // namespace Sierra
