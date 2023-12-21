/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleVofNonConformalSolverAlgorithm.h>
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
// AssembleVofNonConformalSolverAlgorithm - lhs for NC bc (DG)
//                                          used for both edge
//                                          and element
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleVofNonConformalSolverAlgorithm::AssembleVofNonConformalSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  ScalarFieldType *vof)
  : SolverAlgorithm(realm, part, eqSystem),
    vof_(vof),
    coordinates_(NULL),
    exposedAreaVec_(NULL),
    ncVolumeFlowRate_(NULL),
    eta_(realm_.get_nc_alg_upwind_advection() ? 1.0 : 0.0),
    useCurrentNormal_(realm_.get_nc_alg_current_normal())
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<double>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  exposedAreaVec_ = meta_data.get_field<double>(meta_data.side_rank(), "exposed_area_vector");
  ncVolumeFlowRate_ = meta_data.get_field<double>(meta_data.side_rank(), "nc_volume_flow_rate");

  // what do we need ghosted for this alg to work?
  ghostFieldVec_.push_back(&(vof_->field_of_state(stk::mesh::StateNP1)));
  ghostFieldVec_.push_back(coordinates_);
 
  // provide output to user
  if ( useCurrentNormal_ ) 
    NaluEnv::self().naluOutputP0() << "AssembleVofNonConformalSolverAlgorithm::Options: use_current_normal is active" << std::endl;   
  if ( realm_.get_nc_alg_upwind_advection() ) 
    NaluEnv::self().naluOutputP0() << "AssembleVofNonConformalSolverAlgorithm::Options: upwind advective flux is active " << std::endl;   
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleVofNonConformalSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildNonConformalNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleVofNonConformalSolverAlgorithm::execute()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

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

  // mapping for -1:1 -> -0.5:0.5 volume element
  std::vector<double> currentElementIsoParCoords(nDim);
  std::vector<double> opposingElementIsoParCoords(nDim);

  // interpolate nodal values to point-in-elem
  const int sizeOfVofField = 1;
 
  // pointers to fixed values
  double *p_cNx = &cNx[0];
  double *p_oNx = &oNx[0];

  // nodal fields to gather
  std::vector<double> ws_c_face_vof;
  std::vector<double> ws_o_face_vof;
  std::vector<double> ws_c_elem_vof;
  std::vector<double> ws_o_elem_vof;
  std::vector<double> ws_c_elem_coordinates;
  std::vector<double> ws_o_elem_coordinates;
  std::vector<double> ws_o_coordinates;

  // master element data
  std::vector<double> ws_c_det_j;
  std::vector<double> ws_o_det_j;
  std::vector <double > ws_c_general_shape_function;
  std::vector <double > ws_o_general_shape_function;

  // deal with state
  ScalarFieldType &vofNp1 = vof_->field_of_state(stk::mesh::StateNP1);

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
        
        // algorithm related; element
        ws_c_elem_vof.resize(currentNodesPerElement);
        ws_o_elem_vof.resize(opposingNodesPerElement);
        ws_c_elem_coordinates.resize(currentNodesPerElement*nDim);
        ws_o_elem_coordinates.resize(opposingNodesPerElement*nDim);
        ws_c_det_j.resize(1);
        ws_o_det_j.resize(1);
        
        // algorithm related; face
        ws_c_face_vof.resize(currentNodesPerFace);
        ws_o_face_vof.resize(opposingNodesPerFace);
        ws_c_general_shape_function.resize(currentNodesPerFace);
        ws_o_general_shape_function.resize(opposingNodesPerFace);
        ws_o_coordinates.resize(opposingNodesPerFace*nDim);

        // pointers
        double *p_lhs = &lhs[0];
        double *p_rhs = &rhs[0];        

        double *p_c_face_vof = &ws_c_face_vof[0];
        double *p_o_face_vof = &ws_o_face_vof[0];
        double *p_c_elem_vof = &ws_c_elem_vof[0];
        double *p_o_elem_vof = &ws_o_elem_vof[0];
        double *p_c_elem_coordinates = &ws_c_elem_coordinates[0];
        double *p_o_elem_coordinates = &ws_o_elem_coordinates[0];
        double *p_o_coordinates= &ws_o_coordinates[0];

        // me pointers
        double *p_c_general_shape_function = &ws_c_general_shape_function[0];
        double *p_o_general_shape_function = &ws_o_general_shape_function[0];
        
        // populate current face_node_ordinals
        const int *c_face_node_ordinals = meSCSCurrent->side_node_ordinals(currentFaceOrdinal);

        // gather current face data
        stk::mesh::Entity const* current_face_node_rels = bulk_data.begin_nodes(currentFace);
        const int current_num_face_nodes = bulk_data.num_nodes(currentFace);
        for ( int ni = 0; ni < current_num_face_nodes; ++ni ) {
          stk::mesh::Entity node = current_face_node_rels[ni];
          // gather...
          p_c_face_vof[ni] = *stk::mesh::field_data(vofNp1, node);
        }
        
        // populate opposing face_node_ordinals
        const int *o_face_node_ordinals = meSCSOpposing->side_node_ordinals(opposingFaceOrdinal);

        // gather opposing face data
        stk::mesh::Entity const* opposing_face_node_rels = bulk_data.begin_nodes(opposingFace);
        const int opposing_num_face_nodes = bulk_data.num_nodes(opposingFace);
        for ( int ni = 0; ni < opposing_num_face_nodes; ++ni ) {
          stk::mesh::Entity node = opposing_face_node_rels[ni];
          // gather; scalar
          p_o_face_vof[ni] = *stk::mesh::field_data(vofNp1, node);
          // gather; vector
          const double *coords = stk::mesh::field_data(*coordinates_, node);
          for ( int i = 0; i < nDim; ++i ) {
            p_o_coordinates[ni*nDim+i] = coords[i];
          }
        }
        
        // gather current element data; sneak in first of connected nodes
        stk::mesh::Entity const* current_elem_node_rels = bulk_data.begin_nodes(currentElement);
        const int current_num_elem_nodes = bulk_data.num_nodes(currentElement);
        for ( int ni = 0; ni < current_num_elem_nodes; ++ni ) {
          stk::mesh::Entity node = current_elem_node_rels[ni];
          // set connected nodes
          connected_nodes[ni] = node;
          // gather; scalar
          p_c_elem_vof[ni] = *stk::mesh::field_data(vofNp1, node);
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
          p_o_elem_vof[ni] = *stk::mesh::field_data(vofNp1, node);
          // gather; vector
          const double *coords = stk::mesh::field_data(*coordinates_, node);
          const int niNdim = ni*nDim;
          for ( int i = 0; i < nDim; ++i ) {
            p_o_elem_coordinates[niNdim+i] = coords[i];
          }
        }

        // compute opposing normal through master element call, not using oppoing exposed area
        meFCOpposing->general_normal(&opposingIsoParCoords[0], &p_o_coordinates[0], &p_oNx[0]);
        
        // pointer to face data
        const double * c_areaVec = stk::mesh::field_data(*exposedAreaVec_, currentFace);
        const double * ncVolumeFlowRate = stk::mesh::field_data(*ncVolumeFlowRate_, currentFace);
       
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
        
        // interpolate face data; current and opposing...
        double currentVofQBip = 0.0;
        meFCCurrent->interpolatePoint(
          sizeOfVofField,
          &(dgInfo->currentIsoParCoords_[0]),
          &ws_c_face_vof[0],
          &currentVofQBip);
        
        double opposingVofQBip = 0.0;
        meFCOpposing->interpolatePoint(
          sizeOfVofField,
          &(dgInfo->opposingIsoParCoords_[0]),
          &ws_o_face_vof[0],
          &opposingVofQBip);

        // zero lhs/rhs
        for ( int p = 0; p < lhsSize; ++p )
          p_lhs[p] = 0.0;
        for ( int p = 0; p < rhsSize; ++p )
          p_rhs[p] = 0.0;

        // save vdot and |vdot|
        const double tvdot = ncVolumeFlowRate[currentGaussPointId];
        const double abs_tvdot = std::abs(tvdot);
              
        // non conformal advection
        const double ncAdv = tvdot*(currentVofQBip + opposingVofQBip)/2.0 
          + abs_tvdot*(currentVofQBip-opposingVofQBip)/2.0;
       
        // form residual
        const int nn = ipNodeMap[currentGaussPointId];
        p_rhs[nn] -= ncAdv;

        // set-up row for matrix
        const int rowR = nn*totalNodes;
        
        // sensitivities; current face, use general shape function for this single ip
        const double lhsFacC = (eta_*abs_tvdot + tvdot)/2.0;
        meFCCurrent->general_shape_fcn(1, &currentIsoParCoords[0], &ws_c_general_shape_function[0]);
        for ( int ic = 0; ic < currentNodesPerFace; ++ic ) {
          const int icnn = c_face_node_ordinals[ic];
          const double r = p_c_general_shape_function[ic];
          p_lhs[rowR+icnn] += r*lhsFacC;
        }
        
        // sensitivities; opposing face, use general shape function for this single ip
        const double lhsFacO = (eta_*abs_tvdot - tvdot)/2.0;
        meFCOpposing->general_shape_fcn(1, &opposingIsoParCoords[0], &ws_o_general_shape_function[0]);
        for ( int ic = 0; ic < opposingNodesPerFace; ++ic ) {
          const int icnn = o_face_node_ordinals[ic];
          const double r = p_o_general_shape_function[ic];
          p_lhs[rowR+icnn+currentNodesPerElement] -= r*lhsFacO;
        }
        
        apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
