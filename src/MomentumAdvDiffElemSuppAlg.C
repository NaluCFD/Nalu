/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <MomentumAdvDiffElemSuppAlg.h>
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
// MomentumAdvDiffElemSuppAlg - NSO for momentum equation
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MomentumAdvDiffElemSuppAlg::MomentumAdvDiffElemSuppAlg(
  Realm &realm,
  VectorFieldType *velocity,
  ScalarFieldType *viscosity)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    velocityNp1_(NULL),
    coordinates_(NULL),
    viscosity_(viscosity),
    massFlowRate_(NULL),
    nDim_(realm_.meta_data().spatial_dimension()),
    includeDivU_(realm_.get_divU())
{
  // save off fields; for non-BDF2 gather in state N for Nm1 (gamma3_ will be zero)
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  massFlowRate_ = meta_data.get_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "mass_flow_rate_scs");

  // fixed size
  ws_uIp_.resize(nDim_);
}

//--------------------------------------------------------------------------
//-------- elem_resize -----------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumAdvDiffElemSuppAlg::elem_resize(
  MasterElement *meSCS,
  MasterElement */*meSCV*/)
{
  const int nodesPerElement = meSCS->nodesPerElement_;
  const int numScsIp = meSCS->numIntPoints_;

  // resize; geometry
  ws_scs_areav_.resize(numScsIp*nDim_);
  ws_dndx_.resize(nDim_*numScsIp*nodesPerElement);
  ws_deriv_.resize(nDim_*numScsIp*nodesPerElement);
  ws_det_j_.resize(numScsIp);
  ws_shape_function_.resize(numScsIp*nodesPerElement);

  // resize; fields
  ws_uNp1_.resize(nDim_*nodesPerElement);
  ws_coordinates_.resize(nDim_*nodesPerElement);
  ws_viscosity_.resize(nodesPerElement);
  
  // compute shape function
  meSCS->shape_fcn(&ws_shape_function_[0]);
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumAdvDiffElemSuppAlg::setup()
{
  // nothing to extract
}

//--------------------------------------------------------------------------
//-------- elem_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumAdvDiffElemSuppAlg::elem_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity element,
  MasterElement *meSCS,
  MasterElement */*meSCV*/)
{
  // details on this element topo
  const int nodesPerElement = meSCS->nodesPerElement_;
  const int numScsIp = meSCS->numIntPoints_;
  const int *lrscv = meSCS->adjacentNodes();    
  
  // gather
  stk::mesh::Entity const *  node_rels = bulkData_->begin_nodes(element);
  int num_nodes = bulkData_->num_nodes(element);

  // sanity check on num nodes
  ThrowAssert( num_nodes == nodesPerElement );

  for ( int ni = 0; ni < num_nodes; ++ni ) {
    stk::mesh::Entity node = node_rels[ni];
    
    // gather scalars
    ws_viscosity_[ni] = *stk::mesh::field_data(*viscosity_, node);

    // pointers to real data
    const double * uNp1   = stk::mesh::field_data(*velocityNp1_, node );
    const double * coords = stk::mesh::field_data(*coordinates_, node );

    // gather vectors
    const int niNdim = ni*nDim_;
    for ( int i=0; i < nDim_; ++i ) {
      ws_uNp1_[niNdim+i] = uNp1[i];
      ws_coordinates_[niNdim+i] = coords[i];
    }
  }
  
  // compute geometry (AGAIN)...
  double scs_error = 0.0;
  meSCS->determinant(1, &ws_coordinates_[0], &ws_scs_areav_[0], &scs_error);
  
  // compute dndx (AGAIN)...
  meSCS->grad_op(1, &ws_coordinates_[0], &ws_dndx_[0], &ws_deriv_[0], &ws_det_j_[0], &scs_error);

  // ip data for this element
  const double *mdot = stk::mesh::field_data(*massFlowRate_, element);

  for ( int ip = 0; ip < numScsIp; ++ip ) {

    // left and right nodes for this ip
    const int il = lrscv[2*ip];
    const int ir = lrscv[2*ip+1];

    // save off some offsets
    const int ilNdim = il*nDim_;
    const int irNdim = ir*nDim_;
    const int ipNdim = ip*nDim_;
    const int offSetSF = ip*nodesPerElement;

    // save off mdot
    const double tmdot = mdot[ip];

    // compute scs point values; offset to Shape Function; sneak in divU
    double muIp = 0.0;
    double divU = 0.0;
    for ( int i = 0; i < nDim_; ++i )
      ws_uIp_[i] = 0.0;

    for ( int ic = 0; ic < nodesPerElement; ++ic ) {
      const double r = ws_shape_function_[offSetSF+ic];
      muIp += r*ws_viscosity_[ic];
      const int offSetDnDx = nDim_*nodesPerElement*ip + ic*nDim_;
      for ( int j = 0; j < nDim_; ++j ) {
        const double uj = ws_uNp1_[ic*nDim_+j];
        ws_uIp_[j] += r*uj;
        divU += uj*ws_dndx_[offSetDnDx+j];
      }
    }

    // assemble advection; rhs only; add in divU stress (explicit)
    for ( int i = 0; i < nDim_; ++i ) {

      // 2nd order central
      const double uiIp = ws_uIp_[i];

      // total advection; pressure contribution in time term
      const double aflux = tmdot*uiIp;

      // divU stress term
      const double divUstress = 2.0/3.0*muIp*divU*ws_scs_areav_[ipNdim+i]*includeDivU_;

      const int indexL = ilNdim + i;
      const int indexR = irNdim + i;

      // right hand side; L and R
      rhs[indexL] -= aflux + divUstress;
      rhs[indexR] += aflux + divUstress;
    }

    for ( int ic = 0; ic < nodesPerElement; ++ic ) {

      const int icNdim = ic*nDim_;

      // shape function
      const double r = ws_shape_function_[offSetSF+ic];

      // advection and diffusion
      const double lhsfacAdv = r*tmdot;
      
      for ( int i = 0; i < nDim_; ++i ) {

        const int indexL = ilNdim + i;
        const int indexR = irNdim + i;

        const int rowL = indexL*nodesPerElement*nDim_;
        const int rowR = indexR*nodesPerElement*nDim_;

        const int rLiC_i = rowL+icNdim+i;
        const int rRiC_i = rowR+icNdim+i;

        // advection operator lhs; rhs handled above
        // lhs; il then ir
        lhs[rLiC_i] += lhsfacAdv;
        lhs[rRiC_i] -= lhsfacAdv;

        // viscous stress
        const int offSetDnDx = nDim_*nodesPerElement*ip + icNdim;
        double lhs_riC_i = 0.0;
        for ( int j = 0; j < nDim_; ++j ) {

          const double axj = ws_scs_areav_[ipNdim+j];
          const double uj = ws_uNp1_[icNdim+j];

          // -mu*dui/dxj*A_j; fixed i over j loop; see below..
          const double lhsfacDiff_i = -muIp*ws_dndx_[offSetDnDx+j]*axj;
          // lhs; il then ir
          lhs_riC_i += lhsfacDiff_i;

          // -mu*duj/dxi*A_j
          const double lhsfacDiff_j = -muIp*ws_dndx_[offSetDnDx+i]*axj;
          // lhs; il then ir
          lhs[rowL+icNdim+j] += lhsfacDiff_j;
          lhs[rowR+icNdim+j] -= lhsfacDiff_j;
          // rhs; il then ir
          rhs[indexL] -= lhsfacDiff_j*uj;
          rhs[indexR] += lhsfacDiff_j*uj;
        }

        // deal with accumulated lhs and flux for -mu*dui/dxj*Aj
        lhs[rLiC_i] += lhs_riC_i;
        lhs[rRiC_i] -= lhs_riC_i;
        const double ui = ws_uNp1_[icNdim+i];
        rhs[indexL] -= lhs_riC_i*ui;
        rhs[indexR] += lhs_riC_i*ui;
      }
    }
  }
}
  
} // namespace nalu
} // namespace Sierra
