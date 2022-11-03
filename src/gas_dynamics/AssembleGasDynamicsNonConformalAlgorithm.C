/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <gas_dynamics/AssembleGasDynamicsNonConformalAlgorithm.h>
#include <DgInfo.h>
#include <NonConformalInfo.h>
#include <NonConformalManager.h>
#include <Realm.h>
#include <master_element/MasterElement.h>

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

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleGasDynamicsNonConformalAlgorithm - assembles RHS for nc-gd; AUSM+
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleGasDynamicsNonConformalAlgorithm::AssembleGasDynamicsNonConformalAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *density,
  VectorFieldType *momentum,
  VectorFieldType *velocity,
  ScalarFieldType *totalH,
  ScalarFieldType *pressure,
  ScalarFieldType *temperature,
  ScalarFieldType *speedOfSound,
  ScalarFieldType *viscosity,
  ScalarFieldType *thermalCond,
  GenericFieldType *rhsGasDyn)
  : Algorithm(realm, part),
    density_(density),
    momentum_(momentum),
    velocity_(velocity),
    totalH_(totalH),
    pressure_(pressure),
    temperature_(temperature),
    speedOfSound_(speedOfSound),
    viscosity_(viscosity),
    thermalCond_(thermalCond),
    rhsGasDyn_(rhsGasDyn),
    meshVelocity_(NULL),
    exposedAreaVec_(NULL),
    coordinates_(NULL),
    useCurrentNormal_(realm_.get_nc_alg_current_normal()),
    meshMotionFac_(0.0)
{
  // save off mising fields
  stk::mesh::MetaData & metaData = realm_.meta_data();
  if ( realm_.does_mesh_move() ) {
    meshMotionFac_ = 1.0;
    meshVelocity_ = metaData.get_field<double>(stk::topology::NODE_RANK, "mesh_velocity");
  }
  else {
    meshMotionFac_ = 0.0;
    meshVelocity_ = metaData.get_field<double>(stk::topology::NODE_RANK, "velocity");
  }
  exposedAreaVec_ = metaData.get_field<double>(metaData.side_rank(), "exposed_area_vector");
  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // what do we need ghosted for this alg to work?
  ghostFieldVec_.push_back(density_);
  ghostFieldVec_.push_back(momentum_);
  ghostFieldVec_.push_back(velocity_);
  ghostFieldVec_.push_back(totalH_);
  ghostFieldVec_.push_back(pressure_);
  ghostFieldVec_.push_back(temperature_);
  ghostFieldVec_.push_back(speedOfSound_);
  ghostFieldVec_.push_back(viscosity_);
  ghostFieldVec_.push_back(thermalCond_);
  ghostFieldVec_.push_back(coordinates_);
}
  
//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleGasDynamicsNonConformalAlgorithm::execute()
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  stk::mesh::MetaData & metaData = realm_.meta_data();

  // sizes
  const int nDim = metaData.spatial_dimension();
  const int cOffset = nDim;
  const int eOffset = nDim + 1;

  // constants
  const double oneEighth = 1.0/8.0;
  const double threeSixTeenth = 3.0/16.0;

  // ip values; both boundary and opposing surface
  std::vector<double> currentIsoParCoords(nDim);
  std::vector<double> opposingIsoParCoords(nDim);
  std::vector<double> cNx(nDim);
  std::vector<double> oNx(nDim);
  std::vector<double> currentMomentumBip(nDim);
  std::vector<double> opposingMomentumBip(nDim);
  std::vector<double> currentVelocityBip(nDim);
  std::vector<double> opposingVelocityBip(nDim);
  std::vector<double> currentMeshVelocityBip(nDim);
  std::vector<double> currentRhoMeshVelocityBip(nDim);

  // fixed size; gradients
  double ws_dTdxCurrentIp[3];
  double ws_dTdxOpposingIp[3];
  double ws_dudxCurrentIp[3][3];
  double ws_dudxOpposingIp[3][3];
  double ws_tauCurrentIp[3][3];
  double ws_tauOpposingIp[3][3];

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
  std::vector<double> ws_c_density;
  std::vector<double> ws_o_density;
  std::vector<double> ws_c_momentum;
  std::vector<double> ws_o_momentum;
  std::vector<double> ws_c_velocity;
  std::vector<double> ws_o_velocity;
  std::vector<double> ws_c_meshVelocity; // only require current
  std::vector<double> ws_c_totalH;
  std::vector<double> ws_o_totalH;
  std::vector<double> ws_c_pressure;
  std::vector<double> ws_o_pressure;
  std::vector<double> ws_c_temperature;
  std::vector<double> ws_o_temperature;
  std::vector<double> ws_c_speed_of_sound;
  std::vector<double> ws_o_speed_of_sound;
  // diffusion
  std::vector<double> ws_c_viscosity;
  std::vector<double> ws_o_viscosity;
  std::vector<double> ws_c_thermalCond;
  std::vector<double> ws_o_thermalCond;

  std::vector<double> ws_o_coordinates; // only require opposing

  // element; only for diffusion
  std::vector<double> ws_c_elem_velocity;
  std::vector<double> ws_o_elem_velocity;
  std::vector<double> ws_c_elem_temperature;
  std::vector<double> ws_o_elem_temperature;

  // length scale for penalty
  std::vector<double> ws_c_elem_coordinates;
  std::vector<double> ws_o_elem_coordinates;

  // master element data
  std::vector<double> ws_c_dndx;
  std::vector<double> ws_o_dndx;
  std::vector<double> ws_c_det_j;
  std::vector<double> ws_o_det_j;

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
        const int *faceIpNodeMap = meFCCurrent->ipNodeMap();

        // extract some master element info
        const int currentNodesPerFace = meFCCurrent->nodesPerElement_;
        const int opposingNodesPerFace = meFCOpposing->nodesPerElement_;
        const int currentNodesPerElement = meSCSCurrent->nodesPerElement_;
        const int opposingNodesPerElement = meSCSOpposing->nodesPerElement_;
                
        // algorithm related; face
        ws_c_density.resize(currentNodesPerFace);
        ws_o_density.resize(opposingNodesPerFace);
        ws_c_momentum.resize(currentNodesPerFace*nDim);
        ws_o_momentum.resize(opposingNodesPerFace*nDim);
        ws_c_velocity.resize(currentNodesPerFace*nDim);
        ws_o_velocity.resize(opposingNodesPerFace*nDim);
        ws_c_meshVelocity.resize(currentNodesPerFace*nDim);
        ws_c_totalH.resize(currentNodesPerFace);
        ws_o_totalH.resize(opposingNodesPerFace);
        ws_c_pressure.resize(currentNodesPerFace);
        ws_o_pressure.resize(opposingNodesPerFace);
        ws_c_temperature.resize(currentNodesPerFace);
        ws_o_temperature.resize(opposingNodesPerFace);
        ws_c_speed_of_sound.resize(currentNodesPerFace);
        ws_o_speed_of_sound.resize(opposingNodesPerFace);
        ws_c_viscosity.resize(currentNodesPerFace);
        ws_o_viscosity.resize(opposingNodesPerFace);
        ws_c_thermalCond.resize(currentNodesPerFace);
        ws_o_thermalCond.resize(opposingNodesPerFace);
        ws_o_coordinates.resize(opposingNodesPerFace*nDim);

        // algorithm related; element; dndx will be at a single gauss point
        ws_c_elem_velocity.resize(currentNodesPerElement*nDim);
        ws_o_elem_velocity.resize(opposingNodesPerElement*nDim);
        ws_c_elem_temperature.resize(currentNodesPerElement);
        ws_o_elem_temperature.resize(opposingNodesPerElement);
        ws_c_elem_coordinates.resize(currentNodesPerElement*nDim);
        ws_o_elem_coordinates.resize(opposingNodesPerElement*nDim);
        ws_c_dndx.resize(nDim*currentNodesPerElement);
        ws_o_dndx.resize(nDim*opposingNodesPerElement);
        ws_c_det_j.resize(1);
        ws_o_det_j.resize(1);

        // pointers; face
        double *p_c_density = &ws_c_density[0];
        double *p_o_density = &ws_o_density[0];
        double *p_c_momentum = &ws_c_momentum[0];
        double *p_o_momentum = &ws_o_momentum[0];
        double *p_c_velocity = &ws_c_velocity[0];
        double *p_o_velocity= &ws_o_velocity[0];
        double *p_c_meshVelocity = &ws_c_meshVelocity[0];
        double *p_c_totalH = &ws_c_totalH[0];
        double *p_o_totalH = &ws_o_totalH[0];
        double *p_c_pressure = &ws_c_pressure[0];
        double *p_o_pressure = &ws_o_pressure[0];
        double *p_c_temperature = &ws_c_temperature[0];
        double *p_o_temperature = &ws_o_temperature[0];
        double *p_c_speed_of_sound = &ws_c_speed_of_sound[0];
        double *p_o_speed_of_sound = &ws_o_speed_of_sound[0];
        double *p_c_viscosity = &ws_c_viscosity[0];
        double *p_o_viscosity = &ws_o_viscosity[0];
        double *p_c_thermalCond = &ws_c_thermalCond[0];
        double *p_o_thermalCond = &ws_o_thermalCond[0];
        double *p_o_coordinates = &ws_o_coordinates[0];

        // element (for diffusion operator and length scale)
        double *p_c_elem_temperature = &ws_c_elem_temperature[0];
        double *p_o_elem_temperature = &ws_o_elem_temperature[0];
        double *p_c_elem_velocity = &ws_c_elem_velocity[0];
        double *p_o_elem_velocity = &ws_o_elem_velocity[0];
        double *p_c_elem_coordinates = &ws_c_elem_coordinates[0];
        double *p_o_elem_coordinates = &ws_o_elem_coordinates[0];

        // me pointers
        double *p_c_dndx = &ws_c_dndx[0];
        double *p_o_dndx = &ws_o_dndx[0];
        
        // populate current face_node_ordinals
        const int *c_face_node_ordinals = meSCSCurrent->side_node_ordinals(currentFaceOrdinal);

        // gather current face data
        stk::mesh::Entity const* current_face_node_rels = bulkData.begin_nodes(currentFace);
        const int current_num_face_nodes = bulkData.num_nodes(currentFace);
        for ( int ni = 0; ni < current_num_face_nodes; ++ni ) {
          stk::mesh::Entity node = current_face_node_rels[ni];
          // gather; scalar
          p_c_density[ni] = *stk::mesh::field_data(*density_, node);
          p_c_totalH[ni] = *stk::mesh::field_data(*totalH_, node);
          p_c_pressure[ni] = *stk::mesh::field_data(*pressure_, node);
          p_c_temperature[ni] = *stk::mesh::field_data(*temperature_, node);
          p_c_speed_of_sound[ni] = *stk::mesh::field_data(*speedOfSound_, node);
          p_c_viscosity[ni] = *stk::mesh::field_data(*viscosity_, node);
          p_c_thermalCond[ni] = *stk::mesh::field_data(*thermalCond_, node);
          // gather; vector
          const double *momentum = stk::mesh::field_data(*momentum_, node );
          const double *velocity = stk::mesh::field_data(*velocity_, node );
          const double *meshVelocity = stk::mesh::field_data(*meshVelocity_, node );
          for ( int i = 0; i < nDim; ++i ) {
            const int offSet = i*current_num_face_nodes + ni;        
            p_c_momentum[offSet] = momentum[i];
            p_c_velocity[offSet] = velocity[i];
            p_c_meshVelocity[offSet] = meshVelocity[i];
          }
        }
        
        // populate opposing face_node_ordinals
        const int *o_face_node_ordinals = meSCSOpposing->side_node_ordinals(opposingFaceOrdinal);

        // gather opposing face data
        stk::mesh::Entity const* opposing_face_node_rels = bulkData.begin_nodes(opposingFace);
        const int opposing_num_face_nodes = bulkData.num_nodes(opposingFace);
        for ( int ni = 0; ni < opposing_num_face_nodes; ++ni ) {
          stk::mesh::Entity node = opposing_face_node_rels[ni];
          // gather; scalar
          p_o_density[ni] = *stk::mesh::field_data(*density_, node);
          p_o_totalH[ni] = *stk::mesh::field_data(*totalH_, node);
          p_o_pressure[ni] = *stk::mesh::field_data(*pressure_, node);
          p_o_temperature[ni] = *stk::mesh::field_data(*temperature_, node);
          p_o_speed_of_sound[ni] = *stk::mesh::field_data(*speedOfSound_, node);
          p_o_viscosity[ni] = *stk::mesh::field_data(*viscosity_, node);
          p_o_thermalCond[ni] = *stk::mesh::field_data(*thermalCond_, node);
          // gather; vector
          const double *momentum = stk::mesh::field_data(*momentum_, node );
          const double *velocity = stk::mesh::field_data(*velocity_, node );
          const double *coords = stk::mesh::field_data(*coordinates_, node);
          for ( int i = 0; i < nDim; ++i ) {
            const int offSet = i*current_num_face_nodes + ni;        
            p_o_momentum[offSet] = momentum[i];
            p_o_velocity[offSet] = velocity[i];
            p_o_coordinates[ni*nDim+i] = coords[i];
          }
        }

        // gather current element data
        stk::mesh::Entity const* current_elem_node_rels = bulkData.begin_nodes(currentElement);
        const int current_num_elem_nodes = bulkData.num_nodes(currentElement);
        for ( int ni = 0; ni < current_num_elem_nodes; ++ni ) {
          stk::mesh::Entity node = current_elem_node_rels[ni];
          // gather; scalar
          p_c_elem_temperature[ni] = *stk::mesh::field_data(*temperature_, node);
          // gather; vector
          const double *velocity = stk::mesh::field_data(*velocity_, node);
          const double *coords = stk::mesh::field_data(*coordinates_, node);
          const int niNdim = ni*nDim;
          for ( int i = 0; i < nDim; ++i ) {
            p_c_elem_velocity[niNdim+i] = velocity[i];
            p_c_elem_coordinates[niNdim+i] = coords[i];
          }
        }

        // gather opposing element data
        stk::mesh::Entity const* opposing_elem_node_rels = bulkData.begin_nodes(opposingElement);
        const int opposing_num_elem_nodes = bulkData.num_nodes(opposingElement);
        for ( int ni = 0; ni < opposing_num_elem_nodes; ++ni ) {
          stk::mesh::Entity node = opposing_elem_node_rels[ni];
          // gather; scalar
          p_o_elem_temperature[ni] = *stk::mesh::field_data(*temperature_, node);
          // gather; vector
          const double *velocity = stk::mesh::field_data(*velocity_, node);
          const double *coords = stk::mesh::field_data(*coordinates_, node);
          const int niNdim = ni*nDim;
          for ( int i = 0; i < nDim; ++i ) {
            p_o_elem_velocity[niNdim+i] = velocity[i];
            p_o_elem_coordinates[niNdim+i] = coords[i];
          }
        }
        
        // compute opposing normal through master element call, not using oppoing exposed area
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
        
        // current dTdx and dudxIp
        for ( int i = 0; i < nDim; ++i ) {
          ws_dTdxCurrentIp[i] = 0.0;
          for ( int j = 0; j < nDim; ++j ) {
            ws_dudxCurrentIp[i][j] = 0.0;
          }
        }

        for ( int ic = 0; ic < currentNodesPerElement; ++ic ) {
          const int offSetDnDx = ic*nDim; // single intg. point
          const double tempIC = p_c_elem_temperature[ic];
          for ( int i = 0; i < nDim; ++i ) {
            const double dndxj = p_c_dndx[offSetDnDx+i];
            ws_dTdxCurrentIp[i] += dndxj*tempIC;
            for ( int j = 0; j < nDim; ++j ) {
              ws_dudxCurrentIp[i][j] += dndxj*p_c_elem_velocity[ic*nDim+i];
            }
          }
        }

        // opposing dTdx and dudxIp
        for ( int i = 0; i < nDim; ++i ) {
          ws_dTdxOpposingIp[i] = 0.0;
          for ( int j = 0; j < nDim; ++j ) {
            ws_dudxOpposingIp[i][j] = 0.0;
          }
        }

        for ( int ic = 0; ic < opposingNodesPerElement; ++ic ) {
          const int offSetDnDx = ic*nDim; // single intg. point
          const double tempIC = p_o_elem_temperature[ic];
          for ( int i = 0; i < nDim; ++i ) {
            const double dndxj = p_o_dndx[offSetDnDx+i];
            ws_dTdxOpposingIp[i] += dndxj*tempIC;
            for ( int j = 0; j < nDim; ++j ) {
              ws_dudxOpposingIp[i][j] += dndxj*p_o_elem_velocity[ic*nDim+i];
            }
          }
        }

        // divU current and opposing
        double divUc = 0.0;
        double divUo = 0.0;
        for ( int i = 0; i < nDim; ++i ) {
          divUc += ws_dudxCurrentIp[i][i];
          divUo += ws_dudxOpposingIp[i][i];
        }
        
        // interpolate to boundary ips
        double currentRhoBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_density[0],
          &currentRhoBip);
        
        double opposingRhoBip = 0.0;
        meFCOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_density[0],
          &opposingRhoBip);

        double currentTotalHBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_totalH[0],
          &currentTotalHBip);
        
        double opposingTotalHBip = 0.0;
        meFCOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_totalH[0],
          &opposingTotalHBip);

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

        double currentTemperatureBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_temperature[0],
          &currentTemperatureBip);
        
        double opposingTemperatureBip = 0.0;
        meFCOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_temperature[0],
          &opposingTemperatureBip);

        double currentSpeedOfSoundBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_speed_of_sound[0],
          &currentSpeedOfSoundBip);
        
        double opposingSpeedOfSoundBip = 0.0;
        meFCOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_speed_of_sound[0],
          &opposingSpeedOfSoundBip);

        double currentViscosityBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_viscosity[0],
          &currentViscosityBip);
        
        double opposingViscosityBip = 0.0;
        meFCOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_viscosity[0],
          &opposingViscosityBip);

        double currentThermalCondBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_thermalCond[0],
          &currentThermalCondBip);
        
        double opposingThermalCondBip = 0.0;
        meFCOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_thermalCond[0],
          &opposingThermalCondBip);

        // vector quantities
        meFCCurrent->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_momentum[0],
          &currentMomentumBip[0]);
        
        meFCOpposing->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_momentum[0],
          &opposingMomentumBip[0]);

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

        // mesh velocity; only required at current
        meFCCurrent->interpolatePoint(
          sizeOfVectorField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_meshVelocity[0],
          &currentMeshVelocityBip[0]);

        // extract pointers to nearest node fields
        const int nn = faceIpNodeMap[currentGaussPointId];
        stk::mesh::Entity nNode = current_face_node_rels[nn];
        double *rhsGasDyn = stk::mesh::field_data(*rhsGasDyn_, nNode);

        // switch AUSM nomenclature here... Left is current and Right is opposite
        double machNumberC = 0.0;
        double machNumberO = 0.0;
        for ( int j = 0; j < nDim; ++j) {
          const double nj = p_cNx[j];
          const double meanSpeedOfSound = 0.5*(currentSpeedOfSoundBip + opposingSpeedOfSoundBip);
          machNumberC += (currentVelocityBip[j]-meshMotionFac_*currentMeshVelocityBip[j])*nj/meanSpeedOfSound;
          machNumberO += opposingVelocityBip[j]*nj/meanSpeedOfSound;
        }
        
        // AUSM quantities
        const double absMC = std::abs(machNumberC);
        const double absMO = std::abs(machNumberO);
        const double signMC = machNumberC > 0 ? 1.0 : -1.0;
        const double signMO = machNumberO > 0 ? 1.0 : -1.0;
        
        // script{M}+Left and script{M}-Right
        const double scriptMpC = ( absMC >= 1.0) 
          ? 0.5*(machNumberC + absMC)
          : 0.25*std::pow(machNumberC + 1.0, 2.0) 
          + oneEighth*std::pow(machNumberC*machNumberC - 1.0, 2.0);
        
        const double scriptMmO = ( absMO > 1.0) 
          ? 0.5*(machNumberO - absMO)
          : -0.25*std::pow(machNumberO - 1.0, 2.0) 
          - oneEighth*std::pow(machNumberO*machNumberO - 1.0, 2.0);
        
        // script{P}+Left and script{P}-Right
        const double scriptPpC = ( absMC >= 1.0 )
          ? 0.5*(1.0 + signMC) 
          : 0.25*std::pow(machNumberC + 1.0, 2.0)*(2.0 - machNumberC) 
          + threeSixTeenth*machNumberC*std::pow(machNumberC*machNumberC - 1.0, 2.0);
        
        const double scriptPmO = ( absMO >= 1.0 )
          ? 0.5*(1.0 - signMO) 
          : 0.25*std::pow(machNumberO - 1.0, 2.0)*(2.0 + machNumberO) 
          - threeSixTeenth*machNumberO*std::pow(machNumberO*machNumberO - 1.0, 2.0);
        
        // left-right m's and p's
        const double mCO = scriptMpC + scriptMmO;
        const double pCO = scriptPpC*currentPressureBip + scriptPmO*opposingPressureBip;
        
        // momentum first
        for ( int i = 0; i < nDim; ++i ) {
          
          // advective
          const double fmC = currentMomentumBip[i]*currentSpeedOfSoundBip;
          const double fmO = opposingMomentumBip[i]*opposingSpeedOfSoundBip;
          const double fmCO = c_amag*(0.5*mCO*(fmC + fmO) - 0.5*std::abs(mCO)*(fmO - fmC)) 
            + pCO*c_areaVec[currentGaussPointId*nDim+i];
          
          // diffusive
          double dfmNc = 2.0/3.0*currentViscosityBip*divUc*p_cNx[i];
          double dfmNo = 2.0/3.0*opposingViscosityBip*divUo*p_cNx[i];
          ws_tauCurrentIp[i][i] = -2.0/3.0*currentViscosityBip*divUc;
          ws_tauOpposingIp[i][i] = -2.0/3.0*opposingViscosityBip*divUo;
          for ( int j = 0; j < nDim; ++j ) {
            dfmNc += -currentViscosityBip*(ws_dudxCurrentIp[i][j] + ws_dudxCurrentIp[j][i])*p_cNx[j];
            dfmNo += -opposingViscosityBip*(ws_dudxOpposingIp[i][j] + ws_dudxOpposingIp[j][i])*p_oNx[j];
            ws_tauCurrentIp[i][j] += currentViscosityBip*(ws_dudxCurrentIp[i][j] + ws_dudxCurrentIp[j][i]);
            ws_tauOpposingIp[i][j] += opposingViscosityBip*(ws_dudxOpposingIp[i][j] + ws_dudxOpposingIp[j][i]);
          }
          
          // penalty for momentum
          const double penaltyMomIp 
            = (currentViscosityBip*currentInverseLength + opposingViscosityBip*opposingInverseLength)/2.0;

          // advection/diffusion assembly
          rhsGasDyn[i] -= (fmCO + (0.5*(dfmNc - dfmNo) + penaltyMomIp*(currentMomentumBip[i] - opposingMomentumBip[i]))*c_amag);
        }
        
        // continuity 
        const double fcC = currentRhoBip*currentSpeedOfSoundBip;
        const double fcO = opposingRhoBip*opposingSpeedOfSoundBip;
        const double fcCO = c_amag*(0.5*mCO*(fcC + fcO) - 0.5*std::abs(mCO)*(fcO - fcC));
        // no diffusion
        rhsGasDyn[cOffset] -= fcCO;
        
        // total energy advection
        const double feC = currentTotalHBip*currentSpeedOfSoundBip;
        const double feO = opposingTotalHBip*opposingSpeedOfSoundBip;
        const double feCO = c_amag*(0.5*mCO*(feC + feO) - 0.5*std::abs(mCO)*(feO - feC));        

        // total energy diffusion
        double dfeNc = 0.0;
        double dfeNo = 0.0;
        for ( int i = 0; i < nDim; ++i ) {
          dfeNc += -currentThermalCondBip*ws_dTdxCurrentIp[i]*p_cNx[i];
          dfeNo += -opposingThermalCondBip*ws_dTdxOpposingIp[i]*p_oNx[i];
          const double uiIpc = currentVelocityBip[i];
          const double uiIpo = opposingVelocityBip[i];
          for ( int j = 0; j < nDim; ++j ) {
            dfeNc += -uiIpc*ws_tauCurrentIp[i][j]*p_cNx[j];
            dfeNo += -uiIpo*ws_tauOpposingIp[i][j]*p_oNx[j];
          }
        }
      
        // penalty for total energy
        const double penaltyTeIp 
          = (currentThermalCondBip*currentInverseLength + opposingThermalCondBip*opposingInverseLength)/2.0;

        // diffusion, LHS is: d/dxj(qj) - d/dxj(ui*tauij)
        rhsGasDyn[eOffset] -= (feCO + (0.5*(dfeNc - dfeNo) + penaltyTeIp*(currentTemperatureBip - opposingTemperatureBip))*c_amag);
      }     
    }
  }
}


} // namespace nalu
} // namespace Sierra
