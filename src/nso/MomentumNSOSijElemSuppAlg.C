/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <nso/MomentumNSOSijElemSuppAlg.h>
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
// MomentumNSOSijElemSuppAlg - NSO for momentum equation using gijSij form
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
MomentumNSOSijElemSuppAlg<AlgTraits>::MomentumNSOSijElemSuppAlg(
  Realm &realm,
  VectorFieldType *velocity,
  ElemDataRequests& dataPreReqs)
  : SupplementalAlgorithm(realm),
    velocityNp1_(NULL),
    densityNp1_(NULL),
    pressure_(NULL),
    velocityRTM_(NULL),
    coordinates_(NULL),
    Gjp_(NULL),
    lrscv_(realm.get_surface_master_element(AlgTraits::topo_)->adjacentNodes()),
    Cupw_(0.1),
    small_(1.0e-16),
    includeDivU_(realm_.get_divU())
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  velocityNp1_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  densityNp1_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  pressure_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");
 
  // check for mesh motion for proper velocity
  const bool meshMotion = realm_.does_mesh_move();
  if ( meshMotion )
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  Gjp_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx");

  // compute shape function
  MasterElement *meSCS = realm.get_surface_master_element(AlgTraits::topo_);
  meSCS->shape_fcn(&v_shape_function_(0,0));

  // add master elements
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields
  dataPreReqs.add_gathered_nodal_field(*coordinates_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityNp1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityRTM_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*Gjp_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_,1);
  dataPreReqs.add_gathered_nodal_field(*pressure_,1);
  
  // master element data
  dataPreReqs.add_master_element_call(SCS_AREAV);
  dataPreReqs.add_master_element_call(SCS_GRAD_OP);
  dataPreReqs.add_master_element_call(SCS_GIJ);   
}

//--------------------------------------------------------------------------
//-------- element_execute -------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
void
MomentumNSOSijElemSuppAlg<AlgTraits>::element_execute(
  SharedMemView<double **>& lhs,
  SharedMemView<double *>& rhs,
  stk::mesh::Entity element,
  ScratchViews& scratchViews)
{
  SharedMemView<double**>& v_uNp1 = scratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<double**>& v_velocityRTM = scratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<double**>& v_Gjp = scratchViews.get_scratch_view_2D(*Gjp_);
  SharedMemView<double*>& v_rhoNp1 = scratchViews.get_scratch_view_1D(*densityNp1_);
  SharedMemView<double*>& v_pressure = scratchViews.get_scratch_view_1D(*pressure_);

  SharedMemView<double**>& v_scs_areav = scratchViews.scs_areav;
  SharedMemView<double***>& v_dndx = scratchViews.dndx;
  SharedMemView<double***>& v_gijUpper = scratchViews.gijUpper;
  SharedMemView<double***>& v_gijLower = scratchViews.gijLower;

  // compute nodal ke
  for ( int n = 0; n < AlgTraits::nodesPerElement_; ++n ) {
    double ke = 0.0;
    for ( int j = 0; j < AlgTraits::nDim_; ++j )
      ke += v_uNp1(n,j)*v_uNp1(n,j)/2.0;
    v_ke_(n) = ke;
  }
  
  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {

    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];

    // save off some offsets
    const int ilNdim = il*AlgTraits::nDim_;
    const int irNdim = ir*AlgTraits::nDim_;

    // zero out vector that prevail over all components
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      v_rhoVrtmScs_(i) = 0.0;
      v_uNp1Scs_(i) = 0.0;
      v_dpdxScs_(i) = 0.0;
      v_GjpScs_(i) = 0.0;
      v_dkedxScs_(i) = 0.0;
    }
    double rhoScs = 0.0;
    double divU = 0.0;

    // determine scs values of interest
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {

      // save off shape function
      const double r = v_shape_function_(ip,ic);  
      
      // compute scs derivatives and flux derivative
      const double pressureIC = v_pressure(ic);
      const double rhoIC = v_rhoNp1(ic);
      const double keIC = v_ke_(ic);
      rhoScs += r*rhoIC;
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        const double dnj = v_dndx(ip,ic,j);
        const double vrtmj = v_velocityRTM(ic,j);
        const double uNp1 = v_uNp1(ic,j);
        const double Gjp = v_Gjp(ic,j);
     
        v_rhoVrtmScs_(j) += r*rhoIC*vrtmj;
        v_uNp1Scs_(j) += r*uNp1;
        v_dpdxScs_(j) += pressureIC*dnj;
        v_GjpScs_(j) += r*Gjp;
        v_dkedxScs_(j) += keIC*dnj;
        divU += uNp1*dnj;
      }
    }
    
    // form ke residual (based on fine scale momentum residual used in Pstab)
    double keResidual = 0.0;
    for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
      keResidual += v_uNp1Scs_(j)*(v_dpdxScs_(j) - v_GjpScs_(j))/rhoScs/2.0;
    }

    // denominator for nu as well as terms for "upwind" nu
    double gUpperMagGradQ = 0.0;
    double rhoVrtmiGLowerRhoVrtmj = 0.0;
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      const double dkedxScsi = v_dkedxScs_(i);
      const double rhoVrtmi = v_rhoVrtmScs_(i);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        gUpperMagGradQ += dkedxScsi*v_gijUpper(ip,i,j)*v_dkedxScs_(j);
        rhoVrtmiGLowerRhoVrtmj += rhoVrtmi*v_gijLower(ip,i,j)*v_rhoVrtmScs_(j);
      }
    }      
    
    // construct nu from ke residual
    const double nuResidual = rhoScs*std::sqrt((keResidual*keResidual)/(gUpperMagGradQ+small_));
    
    // construct nu from first-order-like approach; SNL-internal write-up (eq 209)
    // for now, only include advection as full set of terms is too diffuse
    const double nuFirstOrder = std::sqrt(rhoVrtmiGLowerRhoVrtmj);
    
    // limit based on first order; Cupw_ is a fudge factor similar to Guermond's approach
    const double nu = std::min(Cupw_*nuFirstOrder, nuResidual);

    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {

      const int icNdim = ic*AlgTraits::nDim_;

      for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
        
        // divU stress term
        const double divUstress = 2.0/3.0*nu*divU*v_scs_areav(ip,i)*v_gijUpper(ip,i,i)*includeDivU_;

        const int indexL = ilNdim + i;
        const int indexR = irNdim + i;
        
        // NSO diffusion-like term; -2*nu*S^*_ij*g^ij = -nu*g^ij*(dui/dxj + duj/dxi - 2/3*nu*rho*divU*del_ij)*Aj
        double lhs_riC_i = 0.0;
        for ( int j = 0; j < AlgTraits::nDim_; ++j ) {

          const double axj = v_scs_areav(ip,j);
          const double uj = v_uNp1(ic,j);

          // -nu*g^ij*dui/dxj*A_j; fixed i over j loop; see below..
          const double lhsfacDiff_i = -nu*v_gijUpper(ip,i,j)*v_dndx(ip,ic,j)*axj;
          // lhs; il then ir
          lhs_riC_i += lhsfacDiff_i;

          // -nu*g^ij*duj/dxi*A_j
          const double lhsfacDiff_j = -nu*v_gijUpper(ip,i,j)*v_dndx(ip,ic,i)*axj;
          // lhs; il then ir
          lhs(indexL,icNdim+j) += lhsfacDiff_j;
          lhs(indexR,icNdim+j) -= lhsfacDiff_j;
          // rhs; il then ir
          rhs(indexL) -= lhsfacDiff_j*uj;
          rhs(indexR) += lhsfacDiff_j*uj;
        }

        // deal with accumulated lhs and flux for -nu*g^ij*dui/dxj*Aj
        lhs(indexL,icNdim+i) += lhs_riC_i;
        lhs(indexR,icNdim+i) -= lhs_riC_i;
        const double ui = v_uNp1(ic,i);
        rhs(indexL) -= lhs_riC_i*ui + divUstress;
        rhs(indexR) += lhs_riC_i*ui + divUstress;
      }
    }
  }
}
  
INSTANTIATE_SUPPLEMENTAL_ALGORITHM(MomentumNSOSijElemSuppAlg);

} // namespace nalu
} // namespace Sierra
