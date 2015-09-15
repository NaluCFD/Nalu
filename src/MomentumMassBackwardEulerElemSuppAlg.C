/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <MomentumMassBackwardEulerElemSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// MomentumMassBackwardEulerElemSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MomentumMassBackwardEulerElemSuppAlg::MomentumMassBackwardEulerElemSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    velocityN_(NULL),
    velocityNp1_(NULL),
    densityN_(NULL),
    densityNp1_(NULL),
    Gjp_(NULL),
    coordinates_(NULL),
    dt_(0.0),
    nDim_(realm_.spatialDimension_),
    useShifted_(false)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  VectorFieldType *velocity = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  velocityN_ = &(velocity->field_of_state(stk::mesh::StateN));
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  Gjp_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // scratch vecs
  uNScv_.resize(nDim_);
  uNp1Scv_.resize(nDim_);
  GjpScv_.resize(nDim_);
}

//--------------------------------------------------------------------------
//-------- elem_resize -----------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumMassBackwardEulerElemSuppAlg::elem_resize(
  MasterElement */*meSCS*/,
  MasterElement *meSCV)
{
  const int nodesPerElement = meSCV->nodesPerElement_;
  const int numScvIp = meSCV->numIntPoints_;

  // resize
  ws_shape_function_.resize(numScvIp*nodesPerElement);
  ws_uN_.resize(nDim_*nodesPerElement);
  ws_uNp1_.resize(nDim_*nodesPerElement);
  ws_Gjp_.resize(nDim_*nodesPerElement);
  ws_rhoN_.resize(nodesPerElement);
  ws_rhoNp1_.resize(nodesPerElement);
  ws_coordinates_.resize(nDim_*nodesPerElement);
  ws_scv_volume_.resize(numScvIp);

  // compute shape function
  if ( useShifted_ )
    meSCV->shifted_shape_fcn(&ws_shape_function_[0]);
  else
    meSCV->shape_fcn(&ws_shape_function_[0]);
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumMassBackwardEulerElemSuppAlg::setup()
{
  dt_ = realm_.get_time_step();
}

//--------------------------------------------------------------------------
//-------- elem_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumMassBackwardEulerElemSuppAlg::elem_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity element,
  MasterElement */*meSCS*/,
  MasterElement *meSCV)
{
  // pointer to ME methods
  const int *ipNodeMap = meSCV->ipNodeMap();
  const int nodesPerElement = meSCV->nodesPerElement_;
  const int numScvIp = meSCV->numIntPoints_;

  // gather
  stk::mesh::Entity const *  node_rels = bulkData_->begin_nodes(element);
  int num_nodes = bulkData_->num_nodes(element);

  // sanity check on num nodes
  ThrowAssert( num_nodes == nodesPerElement );

  for ( int ni = 0; ni < num_nodes; ++ni ) {
    stk::mesh::Entity node = node_rels[ni];
    // pointers to real data
    const double * uN = stk::mesh::field_data(*velocityN_, node );
    const double * uNp1 = stk::mesh::field_data(*velocityNp1_, node );
    const double * Gjp = stk::mesh::field_data(*Gjp_, node );
    const double * coords =  stk::mesh::field_data(*coordinates_, node);
   
    // gather scalars
    ws_rhoN_[ni] = *stk::mesh::field_data(*densityN_, node);
    ws_rhoNp1_[ni] = *stk::mesh::field_data(*densityNp1_, node);

    // gather vectors
    const int niNdim = ni*nDim_;
    for ( int j=0; j < nDim_; ++j ) {
      ws_uN_[niNdim+j] = uN[j];
      ws_uNp1_[niNdim+j] = uNp1[j];
      ws_Gjp_[niNdim+j] = Gjp[j];
      ws_coordinates_[niNdim+j] = coords[j];
    }
  }

  // compute geometry
  double scv_error = 0.0;
  meSCV->determinant(1, &ws_coordinates_[0], &ws_scv_volume_[0], &scv_error);

  for ( int ip = 0; ip < numScvIp; ++ip ) {
      
    // nearest node to ip
    const int nearestNode = ipNodeMap[ip];
    
    // zero out; scalar and vector
    double rhoNScv = 0.0;
    double rhoNp1Scv = 0.0;
    for ( int j =0; j < nDim_; ++j ) {
      uNScv_[j] = 0.0;
      uNp1Scv_[j] = 0.0;
      GjpScv_[j] = 0.0;
    }
      
    const int offSet = ip*nodesPerElement;
    for ( int ic = 0; ic < nodesPerElement; ++ic ) {
      // save off shape function
      const double r = ws_shape_function_[offSet+ic];

      // density
      rhoNScv += r*ws_rhoN_[ic];
      rhoNp1Scv += r*ws_rhoNp1_[ic];

      // velocity
      const int icNdim = ic*nDim_;
      for ( int j = 0; j < nDim_; ++j ) {
        uNScv_[j] += r*ws_uN_[icNdim+j];
        uNp1Scv_[j] += r*ws_uNp1_[icNdim+j];
        GjpScv_[j] += r*ws_Gjp_[icNdim+j];
      }
    }

    // assemble rhs
    const double scV = ws_scv_volume_[ip];
    const int nnNdim = nearestNode*nDim_;
    for ( int i = 0; i < nDim_; ++i ) {
      rhs[nnNdim+i] += 
        -(rhoNp1Scv*uNp1Scv_[i]-rhoNScv*uNScv_[i])*scV/dt_
        -GjpScv_[i]*scV; //- ws_Gjp_[nnNdim+i]*scV;
    }
    
    // manage LHS
    for ( int ic = 0; ic < nodesPerElement; ++ic ) {
      
      const int icNdim = ic*nDim_;
      
      // save off shape function
      const double r = ws_shape_function_[offSet+ic];
      
      const double lhsfac = r*rhoNp1Scv*scV/dt_;
      
      for ( int i = 0; i < nDim_; ++i ) {
        const int indexNN = nnNdim + i;
        const int rowNN = indexNN*nodesPerElement*nDim_;
        const int rNNiC_i = rowNN+icNdim+i;
        lhs[rNNiC_i] += lhsfac;
      }
    }   
  }
}

} // namespace nalu
} // namespace Sierra
