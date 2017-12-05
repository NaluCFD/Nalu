/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleContinuityNonConformalSolverAlgorithm.h>
#include <EquationSystem.h>
#include <DgInfo.h>
#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <NonConformalInfo.h>
#include <NonConformalManager.h>
#include <Realm.h>
#include <TimeIntegrator.h>
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
// AssembleContinuityNonConformalSolverAlgorithm - lhs for NC bc (DG)
//                                                 used for both edge
//                                                 and element; WIP..
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleContinuityNonConformalSolverAlgorithm::AssembleContinuityNonConformalSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  ScalarFieldType *pressure,
  VectorFieldType *Gjp)
  : SolverAlgorithm(realm, part, eqSystem),
    pressure_(pressure),
    Gjp_(Gjp),
    velocity_(NULL),
    meshVelocity_(NULL),
    coordinates_(NULL),
    density_(NULL),
    exposedAreaVec_(NULL),
    meshMotion_(realm_.does_mesh_move()),
    useCurrentNormal_(realm_.get_nc_alg_current_normal()),
    includePstab_(realm_.get_nc_alg_include_pstab() ? 1.0 : 0.0),
    meshFac_(0.0)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  if ( meshMotion_ ) {
    meshVelocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_velocity");
    meshFac_ = 1.0;
  }
  else {   
    meshVelocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
    meshFac_ = 0.0;
  }

  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");  
  
  // what do we need ghosted for this alg to work?
  ghostFieldVec_.push_back(pressure_);
  ghostFieldVec_.push_back(Gjp_);
  ghostFieldVec_.push_back(coordinates_);
  ghostFieldVec_.push_back(velocity_);
  ghostFieldVec_.push_back(density_);

  if ( useCurrentNormal_ )
    NaluEnv::self().naluOutputP0() << "AssembleContinuityNonConformalSolverAlgorithm::Options: use_current_normal is active" << std::endl;
  if ( includePstab_ )
    NaluEnv::self().naluOutputP0() << "AssembleContinuityNonConformalSolverAlgorithm::Options: include_pstab is active" << std::endl;
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleContinuityNonConformalSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildNonConformalNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleContinuityNonConformalSolverAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // continuity equation scales by projection time scale
  const double dt = realm_.get_time_step();
  const double gamma1 = realm_.get_gamma1();
  const double projTimeScale = dt/gamma1;
  
  // deal with interpolation procedure
  const double interpTogether = realm_.get_mdot_interp();
  const double om_interpTogether = 1.0-interpTogether;

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
  std::vector<double> currentRhoVelocityBip(nDim);
  std::vector<double> opposingRhoVelocityBip(nDim);
  std::vector<double> currentRhoMeshVelocityBip(nDim);
  // pressure stabilization
  std::vector<double> currentGjpBip(nDim);
  std::vector<double> opposingGjpBip(nDim);
  std::vector<double> currentDpdxBip(nDim);
  std::vector<double> opposingDpdxBip(nDim);

  // mapping for -1:1 -> -0.5:0.5 volume element
  std::vector<double> currentElementIsoParCoords(nDim);

  // interpolate nodal values to point-in-elem
  const int sizeOfScalarField = 1;
  const int sizeOfVectorField = nDim;
 
  // pointers to fixed values
  double *p_cNx = &cNx[0];
  double *p_oNx = &oNx[0];

  // nodal fields to gather; face
  std::vector<double> ws_c_pressure;
  std::vector<double> ws_c_Gjp;
  std::vector<double> ws_c_velocity;
  std::vector<double> ws_c_mesh_velocity;
  std::vector<double> ws_c_density;
 
  // element
  std::vector<double> ws_o_elem_Gjp;
  std::vector<double> ws_o_elem_velocity;
  std::vector<double> ws_o_elem_density;

  std::vector<double> ws_c_elem_pressure;
  std::vector<double> ws_o_elem_pressure;
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
        stk::mesh::Entity currentElement = dgInfo->currentElement_;
        stk::mesh::Entity opposingElement = dgInfo->opposingElement_;
        const int currentFaceOrdinal = dgInfo->currentFaceOrdinal_;

        // master element; face and volume
        MasterElement * meFCCurrent = dgInfo->meFCCurrent_; 
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
        ws_c_Gjp.resize(currentNodesPerFace*nDim);
        ws_c_velocity.resize(currentNodesPerFace*nDim);
        ws_c_mesh_velocity.resize(currentNodesPerFace*nDim);
        ws_c_density.resize(currentNodesPerFace);
        ws_c_general_shape_function.resize(currentNodesPerFace);
        
        // algorithm related; element; dndx will be at a single gauss point
        ws_o_elem_velocity.resize(opposingNodesPerElement*nDim);
        ws_o_elem_Gjp.resize(opposingNodesPerElement*nDim);
        ws_o_elem_density.resize(opposingNodesPerElement);
        ws_c_elem_pressure.resize(currentNodesPerElement);
        ws_o_elem_pressure.resize(opposingNodesPerElement);
        ws_c_elem_coordinates.resize(currentNodesPerElement*nDim);
        ws_o_elem_coordinates.resize(opposingNodesPerElement*nDim);
        ws_c_dndx.resize(nDim*currentNodesPerElement);
        ws_o_dndx.resize(nDim*opposingNodesPerElement);
        ws_c_det_j.resize(1);
        ws_o_det_j.resize(1);
        ws_o_general_shape_function.resize(opposingNodesPerElement);

        // pointers
        double *p_lhs = &lhs[0];
        double *p_rhs = &rhs[0];
        
        // face
        double *p_c_pressure = &ws_c_pressure[0];
        double *p_c_Gjp = &ws_c_Gjp[0];
        double *p_c_velocity = &ws_c_velocity[0];
        double *p_c_mesh_velocity = &ws_c_mesh_velocity[0];
        double *p_c_density = &ws_c_density[0];

        // element
        double *p_o_elem_Gjp = &ws_o_elem_Gjp[0];
        double *p_o_elem_velocity = &ws_o_elem_velocity[0];
        double *p_o_elem_density = &ws_o_elem_density[0];
        double *p_c_elem_pressure = &ws_c_elem_pressure[0];
        double *p_o_elem_pressure = &ws_o_elem_pressure[0];
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
          // gather; vector
          const double *velocity = stk::mesh::field_data(*velocity_, node );
          const double *meshVelocity = stk::mesh::field_data(*meshVelocity_, node );
          const double *Gjp = stk::mesh::field_data(*Gjp_, node );
          for ( int i = 0; i < nDim; ++i ) {
            const int offSet = i*current_num_face_nodes + ni; 
            p_c_velocity[offSet] = velocity[i];
            p_c_mesh_velocity[offSet] = meshVelocity[i];
            p_c_Gjp[offSet] = Gjp[i];
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
          p_o_elem_density[ni] = *stk::mesh::field_data(*density_,node);
          // gather; vector
          const double *coords = stk::mesh::field_data(*coordinates_, node);
          const double *velocity = stk::mesh::field_data(*velocity_, node );
          const double *Gjp = stk::mesh::field_data(*Gjp_, node );
          const int niNdim = ni*nDim;
          for ( int i = 0; i < nDim; ++i ) {
            p_o_elem_coordinates[niNdim+i] = coords[i];
            const int offSet = i*opposing_num_elem_nodes + ni;        
            p_o_elem_Gjp[offSet] = Gjp[i];
            p_o_elem_velocity[offSet] = velocity[i];
          }
        }
        
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

        // opposing normal is the opposite of the current
        for ( int i = 0; i < nDim; ++i )
          p_oNx[i] = -p_cNx[i];

        // project from side to element; method deals with the -1:1 isInElement range to the proper underlying CVFEM range
        meSCSCurrent->sidePcoords_to_elemPcoords(currentFaceOrdinal, 1, &currentIsoParCoords[0], &currentElementIsoParCoords[0]);
        
        // compute dndx
        double scs_error = 0.0;
        meSCSCurrent->general_face_grad_op(currentFaceOrdinal, &currentElementIsoParCoords[0], 
                                           &p_c_elem_coordinates[0], &p_c_dndx[0], &ws_c_det_j[0], &scs_error);
        meSCSOpposing->general_grad_op(&opposingIsoParCoords[0], 
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

        // opposing inverse length scale; for now, use the same length scale... probably looks like ni*gij*nj
        double opposingInverseLength = currentInverseLength;

        // projected nodal gradient; zero out
        for ( int j = 0; j < nDim; ++j ) {
          currentDpdxBip[j] = 0.0;
          opposingDpdxBip[j] = 0.0;
        }

        // current pressure gradient
        for ( int ic = 0; ic < currentNodesPerElement; ++ic ) {
          const int offSetDnDx = ic*nDim; // single intg. point
          const double pNp1 = p_c_elem_pressure[ic];
          for ( int j = 0; j < nDim; ++j ) {
            const double dndxj = p_c_dndx[offSetDnDx+j];
            currentDpdxBip[j] += dndxj*pNp1;
          }
        }

        // opposing pressure gradient
        for ( int ic = 0; ic < opposingNodesPerElement; ++ic ) {
          const int offSetDnDx = ic*nDim; // single intg. point
          const double pNp1 = p_o_elem_pressure[ic];
          for ( int j = 0; j < nDim; ++j ) {
            const double dndxj = p_o_dndx[offSetDnDx+j];
            opposingDpdxBip[j] += dndxj*pNp1;
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
        meSCSOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_elem_pressure[0],
          &opposingPressureBip);

        // velocity
        meFCCurrent->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_velocity[0],
          &currentVelocityBip[0]);
        
        meSCSOpposing->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_elem_velocity[0],
          &opposingVelocityBip[0]);

        // mesh velocity
        meFCCurrent->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_mesh_velocity[0],
          &currentMeshVelocityBip[0]);

        // projected nodal gradient
        meFCCurrent->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_Gjp[0],
          &currentGjpBip[0]);
        
        meSCSOpposing->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_elem_Gjp[0],
          &opposingGjpBip[0]);

        // density
        double currentDensityBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_density[0],
          &currentDensityBip);
        
        double opposingDensityBip = 0.0;
        meSCSOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_elem_density[0],
          &opposingDensityBip);

        // product of density and velocity; current 
        for ( int ni = 0; ni < current_num_face_nodes; ++ni ) {
          const double density = p_c_density[ni];
          for ( int i = 0; i < nDim; ++i ) {
            const int offSet = i*current_num_face_nodes + ni;        
            p_c_velocity[offSet] *= density;
            p_c_mesh_velocity[offSet] *= density;
          }
        }

        for ( int ni = 0; ni < opposing_num_elem_nodes; ++ni ) {
          const double density = p_o_elem_density[ni];
          for ( int i = 0; i < nDim; ++i ) {
            const int offSet = i*opposing_num_elem_nodes + ni;        
            p_o_elem_velocity[offSet] *= density;
          }
        }

        // interpolate velocity with density scaling
        meFCCurrent->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_velocity[0],
          &currentRhoVelocityBip[0]);
        
        meSCSOpposing->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_elem_velocity[0],
          &opposingRhoVelocityBip[0]);

        // interpolate mesh velocity with density scaling
        meFCCurrent->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_mesh_velocity[0],
          &currentRhoMeshVelocityBip[0]);

        // zero lhs/rhs
        for ( int p = 0; p < lhsSize; ++p )
          p_lhs[p] = 0.0;
        for ( int p = 0; p < rhsSize; ++p )
          p_rhs[p] = 0.0;
                
        const double penaltyIp = projTimeScale*0.5*(currentInverseLength + opposingInverseLength);

        double ncFlux = 0.0;
        double ncMeshFlux = 0.0;
        double ncPstabFlux = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double cRhoVelocity = interpTogether*currentRhoVelocityBip[j] + om_interpTogether*currentDensityBip*currentVelocityBip[j];
          const double oRhoVelocity = interpTogether*opposingRhoVelocityBip[j] + om_interpTogether*opposingDensityBip*opposingVelocityBip[j];
          const double cRhoMeshVelocity = interpTogether*currentRhoMeshVelocityBip[j] + om_interpTogether*currentDensityBip*currentMeshVelocityBip[j];
          ncFlux += 0.5*(cRhoVelocity*p_cNx[j] - oRhoVelocity*p_oNx[j]);
          ncMeshFlux += cRhoMeshVelocity*p_cNx[j];
          const double cPstab = currentDpdxBip[j] - currentGjpBip[j];
          const double oPstab = opposingDpdxBip[j] - opposingGjpBip[j];
          ncPstabFlux += 0.5*(cPstab*p_cNx[j] - oPstab*p_oNx[j]);
        }

        const double mdot = (ncFlux - ncMeshFlux*meshFac_ -includePstab_*projTimeScale*ncPstabFlux 
                             + penaltyIp*(currentPressureBip - opposingPressureBip))*c_amag;
        
        // form residual
        const int nn = ipNodeMap[currentGaussPointId];
        p_rhs[nn] -= mdot/projTimeScale;

        // set-up row for matrix
        const int rowR = nn*totalNodes;
        double lhsFac = penaltyIp*c_amag/projTimeScale;
        
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
          p_lhs[rowR+ic] += 0.5*lhscd*c_amag*includePstab_;
        }

        // sensitivities; opposing element (penalty); use general shape function for this single ip
        meSCSOpposing->general_shape_fcn(1, &opposingIsoParCoords[0], &ws_o_general_shape_function[0]);
        for ( int ic = 0; ic < opposingNodesPerElement; ++ic ) {
          const double r = p_o_general_shape_function[ic];
          p_lhs[rowR+ic+currentNodesPerElement] -= r*lhsFac;
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
          p_lhs[rowR+ic+currentNodesPerElement] -= 0.5*lhscd*c_amag*includePstab_;
        }

        apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
