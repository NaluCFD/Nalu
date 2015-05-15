/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <ComputeMdotNonConformalAlgorithm.h>
#include <Algorithm.h>
#include <DgInfo.h>
#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <NonConformalInfo.h>
#include <NonConformalManager.h>
#include <Realm.h>
#include <SolutionOptions.h>
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
// ComputeMdotNonConformalAlgorithm - compute mdot; both edge and elem
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ComputeMdotNonConformalAlgorithm::ComputeMdotNonConformalAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *pressure,
  ScalarFieldType *ncNormalFlux,
  ScalarFieldType *ncPenalty)
  : Algorithm(realm, part),
    pressure_(pressure),
    ncNormalFlux_(ncNormalFlux),
    ncPenalty_(ncPenalty)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
  ncMassFlowRate_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "nc_mass_flow_rate");
 
  // what do we need ghosted for this alg to work?
  ghostFieldVec_.push_back(pressure_);
  ghostFieldVec_.push_back(ncNormalFlux_);
  ghostFieldVec_.push_back(ncPenalty_);

}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ComputeMdotNonConformalAlgorithm::~ComputeMdotNonConformalAlgorithm()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeMdotNonConformalAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
 
  // ip values; both boundary and opposing surface
  std::vector<double> currentIsoParCoords(nDim-1);
  std::vector<double> opposingIsoParCoords(nDim-1);
  std::vector<double> cNx(nDim);
  std::vector<double> oNx(nDim);

  // interpolate nodal values to point-in-elem
  const int sizeOfScalarField = 1;
 
  // pointers to fixed values
  double *p_cNx = &cNx[0];
  double *p_oNx = &oNx[0];

  // nodal fields to gather
  std::vector<double> ws_c_pressure;
  std::vector<double> ws_o_pressure;
  std::vector<double> ws_c_ncNormalFlux;
  std::vector<double> ws_o_ncNormalFlux;
  std::vector<double> ws_c_ncPenalty;
  std::vector<double> ws_o_ncPenalty;

  std::vector <double > ws_c_general_shape_function;
  std::vector <double > ws_o_general_shape_function;

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
        
        // master element
        MasterElement * meFCCurrent = dgInfo->meFCCurrent_; 
        MasterElement * meFCOpposing = dgInfo->meFCOpposing_;
        
        // local ip, ordinals, etc
        const int currentGaussPointId = dgInfo->currentGaussPointId_;
        currentIsoParCoords = dgInfo->currentIsoParCoords_;
        opposingIsoParCoords = dgInfo->opposingIsoParCoords_;

        // pointer to mdot
        double * ncMassFlowRate = stk::mesh::field_data(*ncMassFlowRate_, currentFace);
        
        // extract some master element info
        const int currentNodesPerFace = meFCCurrent->nodesPerElement_;
        const int opposingNodesPerFace = meFCOpposing->nodesPerElement_;
        
        // algorithm related; face
        ws_c_pressure.resize(currentNodesPerFace);
        ws_o_pressure.resize(opposingNodesPerFace);
        ws_c_ncNormalFlux.resize(currentNodesPerFace);
        ws_o_ncNormalFlux.resize(opposingNodesPerFace);
        ws_c_ncPenalty.resize(currentNodesPerFace);
        ws_o_ncPenalty.resize(opposingNodesPerFace);
        ws_c_general_shape_function.resize(currentNodesPerFace);
        ws_o_general_shape_function.resize(opposingNodesPerFace);
        
        double *p_c_pressure = &ws_c_pressure[0];
        double *p_o_pressure = &ws_o_pressure[0];
        double *p_c_ncNormalFlux = &ws_c_ncNormalFlux[0];
        double *p_o_ncNormalFlux = &ws_o_ncNormalFlux[0];
        double *p_c_ncPenalty = &ws_c_ncPenalty[0];
        double *p_o_ncPenalty = &ws_o_ncPenalty[0];
        
        // general shape function
        double *p_c_general_shape_function = &ws_c_general_shape_function[0];
        double *p_o_general_shape_function = &ws_o_general_shape_function[0];
        
        // gather current face data
        stk::mesh::Entity const* current_face_node_rels = bulk_data.begin_nodes(currentFace);
        const int current_num_face_nodes = bulk_data.num_nodes(currentFace);
        for ( int ni = 0; ni < current_num_face_nodes; ++ni ) {
          stk::mesh::Entity node = current_face_node_rels[ni];
          // gather...
          p_c_pressure[ni] = *stk::mesh::field_data(*pressure_, node);
          p_c_ncNormalFlux[ni] = *stk::mesh::field_data(*ncNormalFlux_, node);
          p_c_ncPenalty[ni] = *stk::mesh::field_data(*ncPenalty_, node);
        }
        
        // gather opposing face data
        stk::mesh::Entity const* opposing_face_node_rels = bulk_data.begin_nodes(opposingFace);
        const int opposing_num_face_nodes = bulk_data.num_nodes(opposingFace);
        for ( int ni = 0; ni < opposing_num_face_nodes; ++ni ) {
          stk::mesh::Entity node = opposing_face_node_rels[ni];
          // gather...
          p_o_pressure[ni] = *stk::mesh::field_data(*pressure_, node);
          p_o_ncNormalFlux[ni] = *stk::mesh::field_data(*ncNormalFlux_, node);
          p_o_ncPenalty[ni] = *stk::mesh::field_data(*ncPenalty_, node);
        }
        
        // pointer to face data
        const double * c_areaVec = stk::mesh::field_data(*exposedAreaVec_, currentFace);
        const double * o_areaVec = stk::mesh::field_data(*exposedAreaVec_, opposingFace);
        
        double c_amag = 0.0;
        double o_amag = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double c_axj = c_areaVec[currentGaussPointId*nDim+j];
          c_amag += c_axj*c_axj;
          // FIXME: choose first area vector on opposing surface? probably need something better for HO
          const double o_axj = o_areaVec[0*nDim+j];
          o_amag += o_axj*o_axj;
        }
        c_amag = std::sqrt(c_amag);
        o_amag = std::sqrt(o_amag);
        
        // now compute normal
        for ( int i = 0; i < nDim; ++i ) {
          p_cNx[i] = c_areaVec[currentGaussPointId*nDim+i]/c_amag;
          p_oNx[i] = o_areaVec[0*nDim+i]/o_amag;  
        }
        
        // interpolate face data
        double currentLambdaBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_ncPenalty[0],
          &currentLambdaBip);
        
        double opposingLambdaBip = 0.0;
        meFCOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_ncPenalty[0],
          &opposingLambdaBip);
        
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

        double currentNcNormalFluxQBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_ncNormalFlux[0],
          &currentNcNormalFluxQBip);

        double opposingNcNormalFluxQBip = 0.0;
        meFCOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_ncNormalFlux[0],
          &opposingNcNormalFluxQBip);
                
        const double penaltyIp = 0.5*(currentLambdaBip + opposingLambdaBip);
        const double ncFlux =  0.5*(currentNcNormalFluxQBip - opposingNcNormalFluxQBip);

        // scatter it
        ncMassFlowRate[currentGaussPointId] = (ncFlux + penaltyIp*(currentPressureBip -opposingPressureBip))*c_amag;

      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
