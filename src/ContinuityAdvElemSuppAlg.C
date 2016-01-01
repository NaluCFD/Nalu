/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ContinuityAdvElemSuppAlg.h>
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
// ContinuityAdvElemSuppAlg - CMM (BDF2) for continuity equation ()p-dof)
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ContinuityAdvElemSuppAlg::ContinuityAdvElemSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    velocityRTM_(NULL),
    Gpdx_(NULL),
    pressure_(NULL),
    densityNp1_(NULL),
    coordinates_(NULL),
    projTimeScale_(1.0),
    nDim_(realm_.spatialDimension_),
    meshMotion_(realm_.does_mesh_move()),
    shiftMdot_(realm_.get_cvfem_shifted_mdot()),
    shiftPoisson_(realm_.get_cvfem_shifted_poisson()),
    reducedSensitivities_(realm_.get_cvfem_reduced_sens_poisson()),
    interpTogether_(realm_.get_mdot_interp()),
    om_interpTogether_(1.0-interpTogether_)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  if ( meshMotion_ )
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  Gpdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx");
  pressure_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // fixed size
  ws_uIp_.resize(nDim_);
  ws_rho_uIp_.resize(nDim_);
  ws_Gpdx_Ip_.resize(nDim_);
  ws_dpdxIp_.resize(nDim_);

  // Implementation details: code is designed to manage the following
  // When shiftPoisson_ is TRUE, reducedSensitivities_ is enforced to be TRUE
  // However, shiftPoisson_ can be FALSE while reducedSensitivities_ is TRUE
}

//--------------------------------------------------------------------------
//-------- elem_resize -----------------------------------------------------
//--------------------------------------------------------------------------
void
ContinuityAdvElemSuppAlg::elem_resize(
  MasterElement *meSCS,
  MasterElement */*meSCV*/)
{
  const int nodesPerElement = meSCS->nodesPerElement_;
  const int numScsIp = meSCS->numIntPoints_;

  // resize
  ws_velocityRTM_.resize(nDim_*nodesPerElement);
  ws_Gpdx_.resize(nDim_*nodesPerElement);
  ws_pressure_.resize(nodesPerElement);
  ws_densityNp1_.resize(nodesPerElement);
  ws_coordinates_.resize(nDim_*nodesPerElement);
  ws_scs_areav_.resize(nDim_*numScsIp);
  ws_dndx_.resize(nDim_*numScsIp*nodesPerElement);
  ws_dndx_lhs_.resize(nDim_*numScsIp*nodesPerElement);
  ws_deriv_.resize(nDim_*numScsIp*nodesPerElement);
  ws_det_j_.resize(numScsIp);
  ws_shape_function_.resize(numScsIp*nodesPerElement);

  // compute shape function
  if ( shiftMdot_ )
    meSCS->shifted_shape_fcn(&ws_shape_function_[0]);
  else
    meSCS->shape_fcn(&ws_shape_function_[0]);
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
ContinuityAdvElemSuppAlg::setup()
{
  const double dt = realm_.get_time_step();
  const double gamma1 = realm_.get_gamma1();
  projTimeScale_ = dt/gamma1;
}

//--------------------------------------------------------------------------
//-------- elem_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
ContinuityAdvElemSuppAlg::elem_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity element,
  MasterElement *meSCS,
  MasterElement */*meSCV*/)
{
  // pointer to ME methods
  const int nodesPerElement = meSCS->nodesPerElement_;
  const int numScsIp = meSCS->numIntPoints_;
  const int *lrscv = meSCS->adjacentNodes();

  // pointers for tricky LHS management
  double *p_dndx = &ws_dndx_[0];
  double *p_dndx_lhs = shiftPoisson_ ? &ws_dndx_[0] : reducedSensitivities_ ? &ws_dndx_lhs_[0] : &ws_dndx_[0];

  // gather
  stk::mesh::Entity const *  node_rels = bulkData_->begin_nodes(element);
  int num_nodes = bulkData_->num_nodes(element);

  // sanity check on num nodes
  ThrowAssert( num_nodes == nodesPerElement );

  for ( int ni = 0; ni < num_nodes; ++ni ) {
    stk::mesh::Entity node = node_rels[ni];

    // pointers to real data
    const double * Gjp    = stk::mesh::field_data(*Gpdx_, node );
    const double * coords = stk::mesh::field_data(*coordinates_, node );
    const double * vrtm   = stk::mesh::field_data(*velocityRTM_, node );

    // gather scalars
    ws_pressure_[ni] = *stk::mesh::field_data(*pressure_, node );
    ws_densityNp1_[ni]  = *stk::mesh::field_data(*densityNp1_, node );

    // gather vectors
    const int niNdim = ni*nDim_;
    for ( int j=0; j < nDim_; ++j ) {
      ws_velocityRTM_[niNdim+j] = vrtm[j];
      ws_Gpdx_[niNdim+j] = Gjp[j];
      ws_coordinates_[niNdim+j] = coords[j];
    }
  }

  // compute geometry
  double scs_error = 0.0;
  meSCS->determinant(1, &ws_coordinates_[0], &ws_scs_areav_[0], &scs_error);

  // compute dndx for residual
  if ( shiftPoisson_ )
    meSCS->shifted_grad_op(1, &ws_coordinates_[0], &p_dndx[0], &ws_deriv_[0], &ws_det_j_[0], &scs_error);
  else
    meSCS->grad_op(1, &ws_coordinates_[0], &p_dndx[0], &ws_deriv_[0], &ws_det_j_[0], &scs_error);

  // compute dndx for LHS
  if ( !shiftPoisson_ && reducedSensitivities_ )
    meSCS->shifted_grad_op(1, &ws_coordinates_[0], &p_dndx_lhs[0], &ws_deriv_[0], &ws_det_j_[0], &scs_error);

  for ( int ip = 0; ip < numScsIp; ++ip ) {

    // left and right nodes for this ip
    const int il = lrscv[2*ip];
    const int ir = lrscv[2*ip+1];

    // corresponding matrix rows
    int rowL = il*nodesPerElement;
    int rowR = ir*nodesPerElement;

    // setup for ip values; sneak in geometry for possible reduced sens
    for ( int j = 0; j < nDim_; ++j ) {
      ws_uIp_[j] = 0.0;
      ws_rho_uIp_[j] = 0.0;
      ws_Gpdx_Ip_[j] = 0.0;
      ws_dpdxIp_[j] = 0.0;
    }
    double rhoIp = 0.0;

    const int offSet = ip*nodesPerElement;
    for ( int ic = 0; ic < nodesPerElement; ++ic ) {

      const double r = ws_shape_function_[offSet+ic];
      const double nodalPressure = ws_pressure_[ic];
      const double nodalRho = ws_densityNp1_[ic];

      rhoIp += r*nodalRho;

      double lhsfac = 0.0;
      const int offSetDnDx = nDim_*nodesPerElement*ip + ic*nDim_;
      for ( int j = 0; j < nDim_; ++j ) {
        ws_Gpdx_Ip_[j] += r*ws_Gpdx_[nDim_*ic+j];
        ws_uIp_[j] += r*ws_velocityRTM_[nDim_*ic+j];
        ws_rho_uIp_[j] += r*nodalRho*ws_velocityRTM_[nDim_*ic+j];
        ws_dpdxIp_[j] += p_dndx[offSetDnDx+j]*nodalPressure;
        lhsfac += -p_dndx_lhs[offSetDnDx+j]*ws_scs_areav_[ip*nDim_+j];
      }

      // assemble to lhs; left
      lhs[rowL+ic] += lhsfac;

      // assemble to lhs; right
      lhs[rowR+ic] -= lhsfac;
    }

    // assemble mdot
    double mdot = 0.0;
    for ( int j = 0; j < nDim_; ++j ) {
      mdot += (interpTogether_*ws_rho_uIp_[j] + om_interpTogether_*rhoIp*ws_uIp_[j]
               - projTimeScale_*(ws_dpdxIp_[j] - ws_Gpdx_Ip_[j]))*ws_scs_areav_[ip*nDim_+j];
    }

    // residual; left and right
    rhs[il] -= mdot/projTimeScale_;
    rhs[ir] += mdot/projTimeScale_;
  }
}
  
} // namespace nalu
} // namespace Sierra
