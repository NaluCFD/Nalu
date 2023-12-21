/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleContinuityVofNonConformalSolverAlgorithm.h>
#include <EquationSystem.h>
#include <DgInfo.h>
#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <NonConformalInfo.h>
#include <NonConformalManager.h>
#include <Realm.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleContinuityVofNonConformalSolverAlgorithm - lhs for NC bc (DG)
//                                                    used for both edge
//                                                    and element
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleContinuityVofNonConformalSolverAlgorithm::AssembleContinuityVofNonConformalSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  ScalarFieldType *pressure,
  VectorFieldType *Gjp,
  const SolutionOptions &solnOpts)
  : SolverAlgorithm(realm, part, eqSystem),
    pressure_(pressure),
    Gjp_(Gjp),
    velocity_(NULL),
    meshVelocity_(NULL),
    coordinates_(NULL),
    density_(NULL),
    vof_(NULL),
    interfaceCurvature_(NULL),
    surfaceTension_(NULL),
    exposedAreaVec_(NULL),
    meshMotion_(realm_.does_mesh_move()),
    useCurrentNormal_(realm_.get_nc_alg_current_normal()),
    includePstab_(realm_.get_nc_alg_include_pstab() ? 1.0 : 0.0),
    meshMotionFac_(0.0),
    n_(solnOpts.localVofN_),
    m_(solnOpts.localVofM_),
    c_(solnOpts.localVofC_),
    buoyancyWeight_(0.0)  
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  velocity_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "velocity");
  if ( meshMotion_ ) {
    meshMotionFac_ = 1.0;
    meshVelocity_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "mesh_velocity");
  }
  else {
    meshMotionFac_ = 0.0;
    meshVelocity_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "velocity");
  }

  coordinates_ = meta_data.get_field<double>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  density_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "density");
  vof_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "volume_of_fluid");
  interfaceCurvature_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "interface_curvature");
  surfaceTension_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "surface_tension");
  exposedAreaVec_ = meta_data.get_field<double>(meta_data.side_rank(), "exposed_area_vector");

  gravity_ = solnOpts.gravity_;

  if( solnOpts.buoyancyPressureStab_ )
    buoyancyWeight_ = 1.0;
  else
    buoyancyWeight_ = 0.0;
  
  // what do we need ghosted for this alg to work?
  ghostFieldVec_.push_back(pressure_);
  ghostFieldVec_.push_back(Gjp_);
  ghostFieldVec_.push_back(coordinates_);
  ghostFieldVec_.push_back(velocity_);
  ghostFieldVec_.push_back(density_);
  ghostFieldVec_.push_back(vof_);
  ghostFieldVec_.push_back(interfaceCurvature_);
  ghostFieldVec_.push_back(surfaceTension_);

  if ( useCurrentNormal_ )
    NaluEnv::self().naluOutputP0() << "AssembleContinuityVofNonConformalSolverAlgorithm::Options: use_current_normal is active" << std::endl;
  if ( includePstab_ )
    NaluEnv::self().naluOutputP0() << "AssembleContinuityVofNonConformalSolverAlgorithm::Options: include_pstab is active" << std::endl;
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleContinuityVofNonConformalSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildNonConformalNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleContinuityVofNonConformalSolverAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // continuity equation scales by projection time scale
  const double dt = realm_.get_time_step();
  const double gamma1 = realm_.get_gamma1();
  const double projTimeScale = dt/gamma1;
  
  // space for LHS/RHS; nodesPerElem*nodesPerElem and nodesPerElem
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<int> scratchIds;
  std::vector<double> scratchVals;
  std::vector<stk::mesh::Entity> connected_nodes;
 
  // ip values; both boundary and opposing surface
  std::vector<double> currentIsoParCoords(nDim);
  std::vector<double> opposingIsoParCoords(nDim);
  std::vector<double> cNx(nDim);
  std::vector<double> oNx(nDim);
  std::vector<double> currentVelocityBip(nDim);
  std::vector<double> opposingVelocityBip(nDim);
  std::vector<double> currentMeshVelocityBip(nDim);

  // pressure stabilization
  std::vector<double> currentGjpBip(nDim);
  std::vector<double> opposingGjpBip(nDim);
  std::vector<double> currentDpdxBip(nDim);
  std::vector<double> opposingDpdxBip(nDim);
  std::vector<double> currentDvofdxBip(nDim);
  std::vector<double> opposingDvofdxBip(nDim);

  // mapping for -1:1 -> -0.5:0.5 volume element
  std::vector<double> currentElementIsoParCoords(nDim);
  std::vector<double> opposingElementIsoParCoords(nDim);

  // interpolate nodal values to point-in-elem
  const int sizeOfScalarField = 1;
  const int sizeOfVectorField = nDim;
 
  // pointers to fixed values
  double *p_cNx = &cNx[0];
  double *p_oNx = &oNx[0];

  // nodal fields to gather; face
  std::vector<double> ws_c_pressure;
  std::vector<double> ws_o_pressure;
  std::vector<double> ws_c_Gjp;
  std::vector<double> ws_o_Gjp;
  std::vector<double> ws_c_velocity;
  std::vector<double> ws_o_velocity;
  std::vector<double> ws_c_meshVelocity; // only require current
  std::vector<double> ws_c_density;
  std::vector<double> ws_o_density;
  std::vector<double> ws_c_vof;
  std::vector<double> ws_o_vof;
  std::vector<double> ws_c_sigma_kappa;
  std::vector<double> ws_o_sigma_kappa;
  std::vector<double> ws_c_surface_tension;
  std::vector<double> ws_o_surface_tension;
  std::vector<double> ws_o_coordinates; // only require opposing

  // element
  std::vector<double> ws_c_elem_pressure;
  std::vector<double> ws_o_elem_pressure;
  std::vector<double> ws_c_elem_vof;
  std::vector<double> ws_o_elem_vof;
  std::vector<double> ws_c_elem_coordinates;
  std::vector<double> ws_o_elem_coordinates;

  // master element data
  std::vector<double> ws_c_dndx;
  std::vector<double> ws_o_dndx;
  std::vector<double> ws_c_det_j;
  std::vector<double> ws_o_det_j;
  std::vector <double > ws_c_general_shape_function;
  std::vector <double > ws_o_general_shape_function;

  // deal with state
  ScalarFieldType &pressureNp1 = pressure_->field_of_state(stk::mesh::StateNP1);

  // parallel communicate ghosted entities
  if ( NULL != realm_.nonConformalManager_->nonConformalGhosting_ )
    stk::mesh::communicate_field_data(*(realm_.nonConformalManager_->nonConformalGhosting_), ghostFieldVec_);

  // iterate nonConformalManager's dgInfoVec
  std::vector<NonConformalInfo *>::iterator ii;
  for( ii=realm_.nonConformalManager_->nonConformalInfoVec_.begin();
       ii!=realm_.nonConformalManager_->nonConformalInfoVec_.end(); ++ii ) {

    // extract vector of DgInfo
    std::vector<std::vector<DgInfo *> > &dgInfoVec = (*ii)->dgInfoVec_;
    
    std::vector<std::vector<DgInfo*> >::iterator idg;
    for( idg=dgInfoVec.begin(); idg!=dgInfoVec.end(); ++idg ) {

      std::vector<DgInfo *> &faceDgInfoVec = (*idg);

      // now loop over all the DgInfo objects on this particular exposed face
      for ( size_t k = 0; k < faceDgInfoVec.size(); ++k ) {

        DgInfo *dgInfo = faceDgInfoVec[k];

        // extract current/opposing face/element
        stk::mesh::Entity currentFace = dgInfo->currentFace_;
        stk::mesh::Entity opposingFace = dgInfo->opposingFace_;
        stk::mesh::Entity currentElement = dgInfo->currentElement_;
        stk::mesh::Entity opposingElement = dgInfo->opposingElement_;
        const int currentFaceOrdinal = dgInfo->currentFaceOrdinal_;
        const int opposingFaceOrdinal = dgInfo->opposingFaceOrdinal_;

        // master element; face and volume
        MasterElement * meFCCurrent = dgInfo->meFCCurrent_; 
        MasterElement * meFCOpposing = dgInfo->meFCOpposing_;
        MasterElement * meSCSCurrent = dgInfo->meSCSCurrent_; 
        MasterElement * meSCSOpposing = dgInfo->meSCSOpposing_;
        
        // local ip, ordinals, etc
        const int currentGaussPointId = dgInfo->currentGaussPointId_;
        currentIsoParCoords = dgInfo->currentIsoParCoords_;
        opposingIsoParCoords = dgInfo->opposingIsoParCoords_;
        
        // mapping from ip to nodes for this ordinal
        const int *ipNodeMap = meSCSCurrent->ipNodeMap(currentFaceOrdinal);

        // extract some master element info
        const int currentNodesPerFace = meFCCurrent->nodesPerElement_;
        const int opposingNodesPerFace = meFCOpposing->nodesPerElement_;
        const int currentNodesPerElement = meSCSCurrent->nodesPerElement_;
        const int opposingNodesPerElement = meSCSOpposing->nodesPerElement_;

        // resize some things; matrix related
        const int totalNodes = currentNodesPerElement + opposingNodesPerElement;
        const int lhsSize = totalNodes*totalNodes;
        const int rhsSize = totalNodes;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connected_nodes.resize(totalNodes);
        
        // algorithm related; face
        ws_c_pressure.resize(currentNodesPerFace);
        ws_o_pressure.resize(opposingNodesPerFace);
        ws_c_Gjp.resize(currentNodesPerFace*nDim);
        ws_o_Gjp.resize(opposingNodesPerFace*nDim);
        ws_c_velocity.resize(currentNodesPerFace*nDim);
        ws_o_velocity.resize(opposingNodesPerFace*nDim);
        ws_c_meshVelocity.resize(currentNodesPerFace*nDim);
        ws_c_density.resize(currentNodesPerFace);
        ws_o_density.resize(opposingNodesPerFace);
        ws_c_vof.resize(currentNodesPerFace);
        ws_o_vof.resize(opposingNodesPerFace);
        ws_c_sigma_kappa.resize(currentNodesPerFace);
        ws_o_sigma_kappa.resize(opposingNodesPerFace);
        ws_o_coordinates.resize(opposingNodesPerFace*nDim);
        ws_c_general_shape_function.resize(currentNodesPerFace);
        ws_o_general_shape_function.resize(opposingNodesPerFace);
        
        // algorithm related; element; dndx will be at a single gauss point
        ws_c_elem_pressure.resize(currentNodesPerElement);
        ws_o_elem_pressure.resize(opposingNodesPerElement);
        ws_c_elem_vof.resize(currentNodesPerElement);
        ws_o_elem_vof.resize(opposingNodesPerElement);
        ws_c_elem_coordinates.resize(currentNodesPerElement*nDim);
        ws_o_elem_coordinates.resize(opposingNodesPerElement*nDim);
        ws_c_dndx.resize(nDim*currentNodesPerElement);
        ws_o_dndx.resize(nDim*opposingNodesPerElement);
        ws_c_det_j.resize(1);
        ws_o_det_j.resize(1);

        // pointers
        double *p_lhs = &lhs[0];
        double *p_rhs = &rhs[0];
        
        // face
        double *p_c_pressure = &ws_c_pressure[0];
        double *p_o_pressure = &ws_o_pressure[0];
        double *p_c_Gjp = &ws_c_Gjp[0];
        double *p_o_Gjp = &ws_o_Gjp[0];
        double *p_c_velocity = &ws_c_velocity[0];
        double *p_o_velocity= &ws_o_velocity[0];
        double *p_c_meshVelocity = &ws_c_meshVelocity[0];
        double *p_c_density = &ws_c_density[0];
        double *p_o_density = &ws_o_density[0];
        double *p_c_vof = &ws_c_vof[0];
        double *p_o_vof = &ws_o_vof[0];
        double *p_c_sigma_kappa = &ws_c_sigma_kappa[0];
        double *p_o_sigma_kappa = &ws_o_sigma_kappa[0];
        double *p_o_coordinates = &ws_o_coordinates[0];

        // element
        double *p_c_elem_pressure = &ws_c_elem_pressure[0];
        double *p_o_elem_pressure = &ws_o_elem_pressure[0];
        double *p_c_elem_vof = &ws_c_elem_vof[0];
        double *p_o_elem_vof = &ws_o_elem_vof[0];
        double *p_c_elem_coordinates = &ws_c_elem_coordinates[0];
        double *p_o_elem_coordinates = &ws_o_elem_coordinates[0];

        // me pointers
        double *p_c_general_shape_function = &ws_c_general_shape_function[0];
        double *p_o_general_shape_function = &ws_o_general_shape_function[0];
        double *p_c_dndx = &ws_c_dndx[0];
        double *p_o_dndx = &ws_o_dndx[0];
        
        // populate current face_node_ordinals
        const int *c_face_node_ordinals = meSCSCurrent->side_node_ordinals(currentFaceOrdinal);

        // gather current face data
        stk::mesh::Entity const* current_face_node_rels = bulk_data.begin_nodes(currentFace);
        const int current_num_face_nodes = bulk_data.num_nodes(currentFace);
        for ( int ni = 0; ni < current_num_face_nodes; ++ni ) {
          stk::mesh::Entity node = current_face_node_rels[ni];          
          // gather; scalar
          p_c_pressure[ni] = *stk::mesh::field_data(pressureNp1, node);
          p_c_density[ni] = *stk::mesh::field_data(*density_, node);
          p_c_vof[ni] = *stk::mesh::field_data(*vof_, node);
          p_c_sigma_kappa[ni] = (*stk::mesh::field_data(*interfaceCurvature_, node))*(*stk::mesh::field_data(*surfaceTension_, node));
          // gather; vector
          const double *velocity = stk::mesh::field_data(*velocity_, node );
          const double *meshVelocity = stk::mesh::field_data(*meshVelocity_, node );
          const double *Gjp = stk::mesh::field_data(*Gjp_, node );
          for ( int i = 0; i < nDim; ++i ) {
            const int offSet = i*current_num_face_nodes + ni; 
            p_c_velocity[offSet] = velocity[i];
            p_c_meshVelocity[offSet] = meshVelocity[i];
            p_c_Gjp[offSet] = Gjp[i];
          }
        }
      
        // populate opposing face_node_ordinals
        const int *o_face_node_ordinals = meSCSOpposing->side_node_ordinals(opposingFaceOrdinal);

        // gather opposing face data
        stk::mesh::Entity const* opposing_face_node_rels = bulk_data.begin_nodes(opposingFace);
        const int opposing_num_face_nodes = bulk_data.num_nodes(opposingFace);
        for ( int ni = 0; ni < opposing_num_face_nodes; ++ni ) {
          stk::mesh::Entity node = opposing_face_node_rels[ni];
          // gather; scalar
          p_o_pressure[ni] = *stk::mesh::field_data(pressureNp1, node);
          p_o_density[ni] = *stk::mesh::field_data(*density_, node);
          p_o_vof[ni] = *stk::mesh::field_data(*vof_, node);
          p_o_sigma_kappa[ni] = (*stk::mesh::field_data(*interfaceCurvature_, node))*(*stk::mesh::field_data(*surfaceTension_, node));
          // gather; vector
          const double *velocity = stk::mesh::field_data(*velocity_, node );
          const double *Gjp = stk::mesh::field_data(*Gjp_, node );
          const double *coords = stk::mesh::field_data(*coordinates_, node);
          for ( int i = 0; i < nDim; ++i ) {
            const int offSet = i*opposing_num_face_nodes + ni;        
            p_o_velocity[offSet] = velocity[i];
            p_o_Gjp[offSet] = Gjp[i];
            p_o_coordinates[ni*nDim+i] = coords[i];
          }
        }
        
        // gather current element data
        stk::mesh::Entity const* current_elem_node_rels = bulk_data.begin_nodes(currentElement);
        const int current_num_elem_nodes = bulk_data.num_nodes(currentElement);
        for ( int ni = 0; ni < current_num_elem_nodes; ++ni ) {
          stk::mesh::Entity node = current_elem_node_rels[ni];          
          // set connected nodes
          connected_nodes[ni] = node;
          // gather; scalar
          p_c_elem_pressure[ni] = *stk::mesh::field_data(pressureNp1, node);
          p_c_elem_vof[ni] = *stk::mesh::field_data(*vof_, node);
          // gather; vector
          const double *coords = stk::mesh::field_data(*coordinates_, node);
          const int niNdim = ni*nDim;
          for ( int i = 0; i < nDim; ++i ) {
            p_c_elem_coordinates[niNdim+i] = coords[i];
          }
        }

        // gather opposing element data; sneak in second connected nodes
        stk::mesh::Entity const* opposing_elem_node_rels = bulk_data.begin_nodes(opposingElement);
        const int opposing_num_elem_nodes = bulk_data.num_nodes(opposingElement);
        for ( int ni = 0; ni < opposing_num_elem_nodes; ++ni ) {
          stk::mesh::Entity node = opposing_elem_node_rels[ni];
          // set connected nodes
          connected_nodes[ni+current_num_elem_nodes] = node;
          // gather; scalar
          p_o_elem_pressure[ni] = *stk::mesh::field_data(pressureNp1, node);
          p_o_elem_vof[ni] = *stk::mesh::field_data(*vof_, node);
          // gather; vector
          const double *coords = stk::mesh::field_data(*coordinates_, node);
          const int niNdim = ni*nDim;
          for ( int i = 0; i < nDim; ++i ) {
            p_o_elem_coordinates[niNdim+i] = coords[i];
          }
        }
        
        // compute opposing normal through master element call, not using opposing exposed area
        meFCOpposing->general_normal(&opposingIsoParCoords[0], &p_o_coordinates[0], &p_oNx[0]);
        
        // pointer to face data
        const double * c_areaVec = stk::mesh::field_data(*exposedAreaVec_, currentFace);
        
        double c_amag = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double c_axj = c_areaVec[currentGaussPointId*nDim+j];
          c_amag += c_axj*c_axj;
        }
        c_amag = std::sqrt(c_amag);

        // now compute normal
        for ( int i = 0; i < nDim; ++i ) {
          p_cNx[i] = c_areaVec[currentGaussPointId*nDim+i]/c_amag;
        }

        // override opposing normal
        if ( useCurrentNormal_ ) {
          for ( int i = 0; i < nDim; ++i )
            p_oNx[i] = -p_cNx[i];
        }

        // project from side to element; method deals with the -1:1 isInElement range to the proper underlying CVFEM range
        meSCSCurrent->sidePcoords_to_elemPcoords(currentFaceOrdinal, 1, &currentIsoParCoords[0], &currentElementIsoParCoords[0]);
        meSCSOpposing->sidePcoords_to_elemPcoords(opposingFaceOrdinal, 1, &opposingIsoParCoords[0], &opposingElementIsoParCoords[0]);
        
        // compute dndx
        double scs_error = 0.0;
        meSCSCurrent->general_face_grad_op(currentFaceOrdinal, &currentElementIsoParCoords[0], 
                                           &p_c_elem_coordinates[0], &p_c_dndx[0], &ws_c_det_j[0], &scs_error);
        meSCSOpposing->general_face_grad_op(opposingFaceOrdinal, &opposingElementIsoParCoords[0], 
                                            &p_o_elem_coordinates[0], &p_o_dndx[0], &ws_o_det_j[0], &scs_error);
        
        // current inverse length scale; can loop over face nodes to avoid "nodesOnFace" array
        double currentInverseLength = 0.0;
        for ( int ic = 0; ic < current_num_face_nodes; ++ic ) {
          const int faceNodeNumber = c_face_node_ordinals[ic];
          const int offSetDnDx = faceNodeNumber*nDim; // single intg. point
          for ( int j = 0; j < nDim; ++j ) {
            const double nxj = p_cNx[j];
            const double dndxj = p_c_dndx[offSetDnDx+j];
            currentInverseLength += dndxj*nxj;
          }
        }

        // opposing inverse length scale; can loop over face nodes to avoid "nodesOnFace" array
        double opposingInverseLength = 0.0;
        for ( int ic = 0; ic < opposing_num_face_nodes; ++ic ) {
          const int faceNodeNumber = o_face_node_ordinals[ic];
          const int offSetDnDx = faceNodeNumber*nDim; // single intg. point
          for ( int j = 0; j < nDim; ++j ) {
            const double nxj = p_oNx[j];
            const double dndxj = p_o_dndx[offSetDnDx+j];
            opposingInverseLength += dndxj*nxj;
          }
        }

        // bip gradients; zero out
        for ( int j = 0; j < nDim; ++j ) {
          currentDpdxBip[j] = 0.0;
          opposingDpdxBip[j] = 0.0;
          currentDvofdxBip[j] = 0.0;
          opposingDvofdxBip[j] = 0.0;
        }

        // current pressure and vof gradient
        for ( int ic = 0; ic < currentNodesPerElement; ++ic ) {
          const int offSetDnDx = ic*nDim; // single intg. point
          const double pNp1 = p_c_elem_pressure[ic];
          const double vofIc = p_c_elem_vof[ic];
          for ( int j = 0; j < nDim; ++j ) {
            const double dndxj = p_c_dndx[offSetDnDx+j];
            currentDpdxBip[j] += dndxj*pNp1;
            currentDvofdxBip[j] += dndxj*vofIc;
          }
        }

        // opposing pressure and vof gradient
        for ( int ic = 0; ic < opposingNodesPerElement; ++ic ) {
          const int offSetDnDx = ic*nDim; // single intg. point
          const double pNp1 = p_o_elem_pressure[ic];
          const double vofIc = p_o_elem_vof[ic];
          for ( int j = 0; j < nDim; ++j ) {
            const double dndxj = p_o_dndx[offSetDnDx+j];
            opposingDpdxBip[j] += dndxj*pNp1;
            opposingDvofdxBip[j] += dndxj*vofIc;
          }
        }

        // interpolate to boundary ips
        double currentPressureBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_pressure[0],
          &currentPressureBip);
        
        double opposingPressureBip = 0.0;
        meFCOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_pressure[0],
          &opposingPressureBip);

        double currentDensityBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_density[0],
          &currentDensityBip);
        
        double opposingDensityBip = 0.0;
        meFCOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_density[0],
          &opposingDensityBip);

        double currentVofBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_vof[0],
          &currentVofBip);
        
        double opposingVofBip = 0.0;
        meFCOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_vof[0],
          &opposingVofBip);

        double currentSigmaKappaBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_sigma_kappa[0],
          &currentSigmaKappaBip);
        
        double opposingSigmaKappaBip = 0.0;
        meFCOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_sigma_kappa[0],
          &opposingSigmaKappaBip);
        
        // projected nodal gradient
        meFCCurrent->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_Gjp[0],
          &currentGjpBip[0]);
        
        meFCOpposing->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_Gjp[0],
          &opposingGjpBip[0]);

        // interpolate velocity
        meFCCurrent->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_velocity[0],
          &currentVelocityBip[0]);
        
        meFCOpposing->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_velocity[0],
          &opposingVelocityBip[0]);

        // interpolate mesh velocity; only current
        meFCCurrent->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_meshVelocity[0],
          &currentMeshVelocityBip[0]);

        // correct for localized approach
        currentSigmaKappaBip *= c_*stk::math::pow(currentVofBip,n_)*stk::math::pow(1.0-currentVofBip,m_);
        opposingSigmaKappaBip *= c_*stk::math::pow(opposingVofBip,n_)*stk::math::pow(1.0-opposingVofBip,m_);

        // zero lhs/rhs
        for ( int p = 0; p < lhsSize; ++p )
          p_lhs[p] = 0.0;
        for ( int p = 0; p < rhsSize; ++p )
          p_rhs[p] = 0.0;
                
        // compute density and penalty factor at the bip
        const double densityBip = 0.5*(currentDensityBip + opposingDensityBip);
        const double penaltyBip = projTimeScale*0.5*(currentInverseLength + opposingInverseLength)/densityBip;

        double ncFlux = 0.0;
        double ncPstabFlux = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double cVelocity = currentVelocityBip[j];
          const double oVelocity = opposingVelocityBip[j];
          const double cMeshVelocity = currentMeshVelocityBip[j];
          ncFlux += 0.5*(cVelocity*p_cNx[j] - oVelocity*p_oNx[j]) - meshMotionFac_*cMeshVelocity*p_cNx[j];
          const double cPstab = (currentDpdxBip[j] 
                                 - buoyancyWeight_*currentDensityBip*gravity_[j] 
                                 + currentSigmaKappaBip*currentDvofdxBip[j])/currentDensityBip - currentGjpBip[j];
          const double oPstab = (opposingDpdxBip[j] 
                                 - buoyancyWeight_*opposingDensityBip*gravity_[j] 
                                 + opposingSigmaKappaBip*opposingDvofdxBip[j])/opposingDensityBip - opposingGjpBip[j];
          ncPstabFlux += 0.5*(cPstab*p_cNx[j] - oPstab*p_oNx[j]);
        }

        const double vdot = (ncFlux - includePstab_*projTimeScale*ncPstabFlux + penaltyBip*(currentPressureBip - opposingPressureBip))*c_amag;
        
        // form residual
        const int nn = ipNodeMap[currentGaussPointId];
        p_rhs[nn] -= vdot/projTimeScale;

        // set-up row for matrix
        const int rowR = nn*totalNodes;
        double lhsFac = penaltyBip*c_amag/projTimeScale;
        
        // sensitivities; current face (penalty); use general shape function for this single ip
        meFCCurrent->general_shape_fcn(1, &currentIsoParCoords[0], &ws_c_general_shape_function[0]);
        for ( int ic = 0; ic < currentNodesPerFace; ++ic ) {
          const int icnn = c_face_node_ordinals[ic];
          const double r = p_c_general_shape_function[ic];
          p_lhs[rowR+icnn] += r*lhsFac;
        }
        
        // sensitivities; current element (diffusion)
        for ( int ic = 0; ic < currentNodesPerElement; ++ic ) {
          const int offSetDnDx = ic*nDim; // single intg. point
          double lhscd = 0.0;
          for ( int j = 0; j < nDim; ++j ) {
            const double nxj = p_cNx[j];
            const double dndxj = p_c_dndx[offSetDnDx+j];
            lhscd -= dndxj*nxj;
          }
          p_lhs[rowR+ic] += 0.5*lhscd*c_amag/currentDensityBip*includePstab_;
        }

        // sensitivities; opposing face (penalty); use general shape function for this single ip
        meFCOpposing->general_shape_fcn(1, &opposingIsoParCoords[0], &ws_o_general_shape_function[0]);
        for ( int ic = 0; ic < opposingNodesPerFace; ++ic ) {
          const int icnn = o_face_node_ordinals[ic];
          const double r = p_o_general_shape_function[ic];
          p_lhs[rowR+icnn+currentNodesPerElement] -= r*lhsFac;
        }
        
        // sensitivities; opposing element (diffusion)
        for ( int ic = 0; ic < opposingNodesPerElement; ++ic ) {
          const int offSetDnDx = ic*nDim; // single intg. point
          double lhscd = 0.0;
          for ( int j = 0; j < nDim; ++j ) {
            const double nxj = p_oNx[j];
            const double dndxj = p_o_dndx[offSetDnDx+j];
            lhscd -= dndxj*nxj;
          }
          p_lhs[rowR+ic+currentNodesPerElement] -= 0.5*lhscd*c_amag/opposingDensityBip*includePstab_;
        }

        apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
