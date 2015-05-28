/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleScalarNonConformalSolverAlgorithm.h>
#include <EquationSystem.h>
#include <DgInfo.h>
#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <NonConformalInfo.h>
#include <NonConformalManager.h>
#include <Realm.h>
#include <SolutionOptions.h>
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
// AssembleScalarNonConformalSolverAlgorithm - lhs for NC bc (DG)
//                                                     used for both edge
//                                                     and element; WIP..
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleScalarNonConformalSolverAlgorithm::AssembleScalarNonConformalSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  ScalarFieldType *scalarQ,
  ScalarFieldType *ncNormalFlux,
  ScalarFieldType *ncPenalty)
  : SolverAlgorithm(realm, part, eqSystem),
    scalarQ_(scalarQ),
    ncNormalFlux_(ncNormalFlux),
    ncPenalty_(ncPenalty),
    exposedAreaVec_(NULL),
    ncMassFlowRate_(NULL),
    robinStyle_(false),
    dsFactor_(1.0)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");  
  ncMassFlowRate_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "nc_mass_flow_rate");

  // what do we need ghosted for this alg to work?
  ghostFieldVec_.push_back(&(scalarQ_->field_of_state(stk::mesh::StateNP1)));
  ghostFieldVec_.push_back(ncNormalFlux_);
  ghostFieldVec_.push_back(ncPenalty_);
 
  // specific algorithm options
  NonConformalAlgType algType = realm_.get_nc_alg_type();
  switch ( algType ) {
    case NC_ALG_TYPE_DG:
      dsFactor_ = 1.0;
      robinStyle_ = false;
      break;
      
    case NC_ALG_TYPE_DS: 
      dsFactor_ = 0.0;
      // robinStyle_ does not matter here..
      break;
    
    case NC_ALG_TYPE_RB:
      dsFactor_ = 1.0;
      robinStyle_ = true;
      
    default:
      // nothing to do... parsing should have caught this...
      break;
  }

  NaluEnv::self().naluOutputP0() << "NC options: dsFactor/robinStyle: " << dsFactor_ << " " << robinStyle_ << std::endl;
  
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarNonConformalSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildNonConformalNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarNonConformalSolverAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // space for LHS/RHS; nodesPerElem*nodesPerElem and nodesPerElem
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<stk::mesh::Entity> connected_nodes;
 
  // ip values; both boundary and opposing surface
  std::vector<double> currentIsoParCoords(nDim-1);
  std::vector<double> opposingIsoParCoords(nDim-1);
  std::vector<double> cNx(nDim);
  std::vector<double> oNx(nDim);

  // interpolate nodal values to point-in-elem
  const int sizeOfScalarField = 1;
  const int sizeOfVectorField = nDim;
 
  // pointers to fixed values
  double *p_cNx = &cNx[0];
  double *p_oNx = &oNx[0];

  // nodal fields to gather
  std::vector<double> ws_c_face_scalarQ;
  std::vector<double> ws_o_face_scalarQ;
  std::vector<double> ws_c_ncNormalFlux;
  std::vector<double> ws_o_ncNormalFlux;
  std::vector<double> ws_c_ncPenalty;
  std::vector<double> ws_o_ncPenalty;

  std::vector <double > ws_c_general_shape_function;
  std::vector <double > ws_o_general_shape_function;

  // deal with state
  ScalarFieldType &scalarQNp1 = scalarQ_->field_of_state(stk::mesh::StateNP1);

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
        
        // extract some master element info
        const int currentNodesPerFace = meFCCurrent->nodesPerElement_;
        const int opposingNodesPerFace = meFCOpposing->nodesPerElement_;
        
        // resize some things; matrix related
        const int lhsSize = (currentNodesPerFace+opposingNodesPerFace)*(currentNodesPerFace+opposingNodesPerFace);
        const int rhsSize = currentNodesPerFace+opposingNodesPerFace;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        connected_nodes.resize(currentNodesPerFace+opposingNodesPerFace);
        
        // algorithm related; element (n/a)
        
        // algorithm related; face
        ws_c_face_scalarQ.resize(currentNodesPerFace);
        ws_o_face_scalarQ.resize(opposingNodesPerFace);
        ws_c_ncNormalFlux.resize(currentNodesPerFace);
        ws_o_ncNormalFlux.resize(opposingNodesPerFace);
        ws_c_ncPenalty.resize(currentNodesPerFace);
        ws_o_ncPenalty.resize(opposingNodesPerFace);
        ws_c_general_shape_function.resize(currentNodesPerFace);
        ws_o_general_shape_function.resize(opposingNodesPerFace);
        
        // pointers
        double *p_lhs = &lhs[0];
        double *p_rhs = &rhs[0];
        
        double *p_c_face_scalarQ = &ws_c_face_scalarQ[0];
        double *p_o_face_scalarQ = &ws_o_face_scalarQ[0];
        double *p_c_ncNormalFlux = &ws_c_ncNormalFlux[0];
        double *p_o_ncNormalFlux = &ws_o_ncNormalFlux[0];
        double *p_c_ncPenalty = &ws_c_ncPenalty[0];
        double *p_o_ncPenalty = &ws_o_ncPenalty[0];
               
        // general shape function
        double *p_c_general_shape_function = &ws_c_general_shape_function[0];
        double *p_o_general_shape_function = &ws_o_general_shape_function[0];
        
        // gather current face data; sneak in first of connected nodes
        stk::mesh::Entity const* current_face_node_rels = bulk_data.begin_nodes(currentFace);
        const int current_num_face_nodes = bulk_data.num_nodes(currentFace);
        for ( int ni = 0; ni < current_num_face_nodes; ++ni ) {
          stk::mesh::Entity node = current_face_node_rels[ni];
          // set connected nodes
          connected_nodes[ni] = node;
          // gather...
          p_c_face_scalarQ[ni] = *stk::mesh::field_data(scalarQNp1, node);
          p_c_ncNormalFlux[ni] = *stk::mesh::field_data(*ncNormalFlux_, node);
          p_c_ncPenalty[ni] = *stk::mesh::field_data(*ncPenalty_, node);
        }
        
        // gather opposing face data; sneak in second of connected nodes
        stk::mesh::Entity const* opposing_face_node_rels = bulk_data.begin_nodes(opposingFace);
        const int opposing_num_face_nodes = bulk_data.num_nodes(opposingFace);
        for ( int ni = 0; ni < opposing_num_face_nodes; ++ni ) {
          stk::mesh::Entity node = opposing_face_node_rels[ni];
          // set connected nodes
          connected_nodes[ni+current_num_face_nodes] = node;
          // gather...
          p_o_face_scalarQ[ni] = *stk::mesh::field_data(scalarQNp1, node);
          p_o_ncNormalFlux[ni] = *stk::mesh::field_data(*ncNormalFlux_, node);
          p_o_ncPenalty[ni] = *stk::mesh::field_data(*ncPenalty_, node);
        }
        
        // pointer to face data
        const double * c_areaVec = stk::mesh::field_data(*exposedAreaVec_, currentFace);
        const double * o_areaVec = stk::mesh::field_data(*exposedAreaVec_, opposingFace);
        const double * ncMassFlowRate = stk::mesh::field_data(*ncMassFlowRate_, currentFace);
       
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
        
        double currentScalarQBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_face_scalarQ[0],
          &currentScalarQBip);
        
        double opposingScalarQBip = 0.0;
        meFCOpposing->interpolatePoint(
          sizeOfScalarField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_face_scalarQ[0],
          &opposingScalarQBip);

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
                
        // zero lhs/rhs
        for ( int p = 0; p < lhsSize; ++p )
          p_lhs[p] = 0.0;
        for ( int p = 0; p < rhsSize; ++p )
          p_rhs[p] = 0.0;

        // compute penalty; did not include mass flow rate
        const double penaltyIp = 0.5*(currentLambdaBip + opposingLambdaBip) + std::abs(ncMassFlowRate[currentGaussPointId])/2.0;
        const double ncDiffFlux =  robinStyle_ ? -opposingNcNormalFluxQBip : 0.5*(currentNcNormalFluxQBip - opposingNcNormalFluxQBip);
        const double ncAdv = robinStyle_ ? ncMassFlowRate[currentGaussPointId]*opposingScalarQBip 
          : 0.5*ncMassFlowRate[currentGaussPointId]*(currentScalarQBip + opposingScalarQBip);
       
        // form residual
        const int nn = currentGaussPointId;
        p_rhs[nn] -= ((dsFactor_*ncDiffFlux + penaltyIp*(currentScalarQBip-opposingScalarQBip))*c_amag + ncAdv);

        // set-up row for matrix
        const int rowR = nn*(currentNodesPerFace+opposingNodesPerFace);
        double lhsFac = penaltyIp*c_amag;
        
        // sensitivities; current face; use general shape function for this single ip
        meFCCurrent->general_shape_fcn(1, &currentIsoParCoords[0], &ws_c_general_shape_function[0]);
        for ( int ic = 0; ic < currentNodesPerFace; ++ic ) {
          const double r = p_c_general_shape_function[ic];
          p_lhs[rowR+ic] += r*lhsFac;
        }
        
        // sensitivities; opposing face; use general shape function for this single ip
        meFCOpposing->general_shape_fcn(1, &opposingIsoParCoords[0], &ws_o_general_shape_function[0]);
        for ( int ic = 0; ic < opposingNodesPerFace; ++ic ) {
          const double r = p_o_general_shape_function[ic];
          p_lhs[rowR+ic+currentNodesPerFace] -= r*lhsFac;
        }
        
        apply_coeff(connected_nodes, rhs, lhs, __FILE__);
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
