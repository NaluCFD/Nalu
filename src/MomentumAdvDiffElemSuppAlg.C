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

// template and scratch space
#include <BuildTemplates.h>
#include <ScratchViews.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>

// topology
#include <stk_topology/topology.hpp>

// Kokkos
#include <Kokkos_Core.hpp>

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
template<class AlgTraits>
MomentumAdvDiffElemSuppAlg<AlgTraits>::MomentumAdvDiffElemSuppAlg(
  Realm &realm,
  VectorFieldType *velocity,
  ScalarFieldType *viscosity,
  ElemDataRequests& dataPreReqs)
  : SupplementalAlgorithm(realm),
    velocityNp1_(NULL),
    coordinates_(NULL),
    viscosity_(viscosity),
    massFlowRate_(NULL),
    lrscv_(realm.get_surface_master_element(AlgTraits::topo_)->adjacentNodes()),
    includeDivU_(realm_.get_divU())
{
  // save off fields; for non-BDF2 gather in state N for Nm1 (gamma3_ will be zero)
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  massFlowRate_ = meta_data.get_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "mass_flow_rate_scs");

  // compute shape function; do we want to push this to dataPreReqs?
  MasterElement *meSCS = realm.get_surface_master_element(AlgTraits::topo_);
  meSCS->shape_fcn(&v_shape_function_(0,0));

  // add master elements
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields and data; mdot not gathered as element data
  dataPreReqs.add_gathered_nodal_field(*coordinates_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocity, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*viscosity_, 1);
  dataPreReqs.add_master_element_call(SCS_AREAV);
  dataPreReqs.add_master_element_call(SCS_GRAD_OP);
}

//--------------------------------------------------------------------------
//-------- element_execute -------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
void
MomentumAdvDiffElemSuppAlg<AlgTraits>::element_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity element,
  ScratchViews& scratchViews)
{
  SharedMemView<double**>& v_uNp1 = scratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<double*>& v_viscosity = scratchViews.get_scratch_view_1D(*viscosity_);

  SharedMemView<double**>& v_scs_areav = scratchViews.scs_areav;
  SharedMemView<double***>& v_dndx = scratchViews.dndx;

  // ip data for this element
  const double *mdot = stk::mesh::field_data(*massFlowRate_, element);

  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {

    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];

    // save off some offsets
    const int ilNdim = il*AlgTraits::nDim_;
    const int irNdim = ir*AlgTraits::nDim_;
  
    // save off mdot
    const double tmdot = mdot[ip];

    // compute scs point values; sneak in divU
    double muIp = 0.0;
    double divU = 0.0;
    for ( int i = 0; i < AlgTraits::nDim_; ++i )
      v_uIp_(i) = 0.0;

    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const double r = v_shape_function_(ip,ic);
      muIp += r*v_viscosity(ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        const double uj = v_uNp1(ic,j);
        v_uIp_[j] += r*uj;
        divU += uj*v_dndx(ip,ic,j);
      }
    }

    // assemble advection; rhs only; add in divU stress (explicit)
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {

      // 2nd order central
      const double uiIp = v_uIp_(i);

      // total advection; (pressure contribution in time term)
      const double aflux = tmdot*uiIp;

      // divU stress term
      const double divUstress = 2.0/3.0*muIp*divU*v_scs_areav(ip,i)*includeDivU_;

      const int indexL = ilNdim + i;
      const int indexR = irNdim + i;

      // right hand side; L and R
      rhs[indexL] -= aflux + divUstress;
      rhs[indexR] += aflux + divUstress;
    }

    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {

      const int icNdim = ic*AlgTraits::nDim_;

      // shape function
      const double r = v_shape_function_(ip,ic);

      // advection and diffusion
      const double lhsfacAdv = r*tmdot;
      
      for ( int i = 0; i < AlgTraits::nDim_; ++i ) {

        const int indexL = ilNdim + i;
        const int indexR = irNdim + i;

        const int rowL = indexL*AlgTraits::nodesPerElement_*AlgTraits::nDim_;
        const int rowR = indexR*AlgTraits::nodesPerElement_*AlgTraits::nDim_;

        const int rLiC_i = rowL+icNdim+i;
        const int rRiC_i = rowR+icNdim+i;

        // advection operator lhs; rhs handled above
        // lhs; il then ir
        lhs[rLiC_i] += lhsfacAdv;
        lhs[rRiC_i] -= lhsfacAdv;

        // viscous stress
        double lhs_riC_i = 0.0;
        for ( int j = 0; j < AlgTraits::nDim_; ++j ) {

          const double axj = v_scs_areav(ip,j);
          const double uj = v_uNp1(ic,j);

          // -mu*dui/dxj*A_j; fixed i over j loop; see below..
          const double lhsfacDiff_i = -muIp*v_dndx(ip,ic,j)*axj;
          // lhs; il then ir
          lhs_riC_i += lhsfacDiff_i;

          // -mu*duj/dxi*A_j
          const double lhsfacDiff_j = -muIp*v_dndx(ip,ic,i)*axj;
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
        const double ui = v_uNp1(ic,i);
        rhs[indexL] -= lhs_riC_i*ui;
        rhs[indexR] += lhs_riC_i*ui;
      }
    }
  }
}
  
INSTANTIATE_SUPPLEMENTAL_ALGORITHM(MomentumAdvDiffElemSuppAlg);

} // namespace nalu
} // namespace Sierra
