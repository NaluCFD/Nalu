/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <nso/MomentumNSOElemSuppAlg.h>
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
// MomentumNSOElemSuppAlg - NSO for momentum equation
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
MomentumNSOElemSuppAlg<AlgTraits>::MomentumNSOElemSuppAlg(
  Realm &realm,
  VectorFieldType *velocity,
  GenericFieldType *Gju,
  ScalarFieldType *viscosity,
  const double fourthFac,
  const double altResFac,
  ElemDataRequests& dataPreReqs)
  : SupplementalAlgorithm(realm),
    velocityNm1_(NULL),
    velocityN_(NULL),
    velocityNp1_(NULL),
    densityNm1_(NULL),
    densityN_(NULL),
    densityNp1_(NULL),
    pressure_(NULL),
    velocityRTM_(NULL),
    coordinates_(NULL),
    viscosity_(viscosity),
    Gju_(Gju),
    lrscv_(realm.get_surface_master_element(AlgTraits::topo_)->adjacentNodes()),
    dt_(0.0),
    gamma1_(0.0),
    gamma2_(0.0),
    gamma3_(0.0),
    Cupw_(0.1),
    small_(1.0e-16),
    fourthFac_(fourthFac),
    altResFac_(altResFac),
    om_altResFac_(1.0-altResFac),
    nonConservedForm_(0.0),
    includeDivU_(realm_.get_divU())
{
  // save off fields; for non-BDF2 gather in state N for Nm1 (gamma3_ will be zero)
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  velocityNm1_ = realm_.number_of_states() == 2 ? &(velocity->field_of_state(stk::mesh::StateN)) : &(velocity->field_of_state(stk::mesh::StateNM1));
  velocityN_ = &(velocity->field_of_state(stk::mesh::StateN));
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNm1_ = realm_.number_of_states() == 2 ? &(density->field_of_state(stk::mesh::StateN)) : &(density->field_of_state(stk::mesh::StateNM1));
  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  pressure_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");
 
  // check for mesh motion for proper velocity
  const bool meshMotion = realm_.does_mesh_move();
  if ( meshMotion )
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // compute shape function; do we want to push this to dataPreReqs?
  MasterElement *meSCS = realm.get_surface_master_element(AlgTraits::topo_);
  meSCS->shape_fcn(&v_shape_function_(0,0));

  // add master elements
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields
  dataPreReqs.add_gathered_nodal_field(*Gju_, AlgTraits::nDim_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*coordinates_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityNm1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityN_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityNp1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityRTM_, AlgTraits::nDim_);
  
  dataPreReqs.add_gathered_nodal_field(*densityNm1_,1);
  dataPreReqs.add_gathered_nodal_field(*densityN_,1);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_,1);
  dataPreReqs.add_gathered_nodal_field(*viscosity_,1);
  dataPreReqs.add_gathered_nodal_field(*pressure_,1);
  
  // master element data
  dataPreReqs.add_master_element_call(SCS_AREAV);
  dataPreReqs.add_master_element_call(SCS_GRAD_OP);
  dataPreReqs.add_master_element_call(SCS_GIJ); 

  // initialize kd
  for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
    for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
      v_kd_(i,j) = (i == j) ? 1.0 : 0.0;
    }
  }
}
  
//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
void
MomentumNSOElemSuppAlg<AlgTraits>::setup()
{
  dt_ = realm_.get_time_step();
  gamma1_ = realm_.get_gamma1();
  gamma2_ = realm_.get_gamma2();
  gamma3_ = realm_.get_gamma3();
}

//--------------------------------------------------------------------------
//-------- element_execute -------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
void
MomentumNSOElemSuppAlg<AlgTraits>::element_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity element,
  ScratchViews& scratchViews)
{
  SharedMemView<double***>& v_Gju = scratchViews.get_scratch_view_3D(*Gju_);
  SharedMemView<double**>& v_uNm1 = scratchViews.get_scratch_view_2D(*velocityNm1_);
  SharedMemView<double**>& v_uN = scratchViews.get_scratch_view_2D(*velocityN_);
  SharedMemView<double**>& v_uNp1 = scratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<double**>& v_velocityRTM = scratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<double*>& v_rhoNm1 = scratchViews.get_scratch_view_1D(*densityNm1_);
  SharedMemView<double*>& v_rhoN = scratchViews.get_scratch_view_1D(*densityN_);
  SharedMemView<double*>& v_rhoNp1 = scratchViews.get_scratch_view_1D(*densityNp1_);
  SharedMemView<double*>& v_viscosity = scratchViews.get_scratch_view_1D(*viscosity_);
  SharedMemView<double*>& v_pressure = scratchViews.get_scratch_view_1D(*pressure_);

  SharedMemView<double**>& v_scs_areav = scratchViews.scs_areav;
  SharedMemView<double***>& v_dndx = scratchViews.dndx;
  SharedMemView<double***>& v_gijUpper = scratchViews.gijUpper;
  SharedMemView<double***>& v_gijLower = scratchViews.gijLower;

  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {

    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];

    // save off some offsets
    const int ilNdim = il*AlgTraits::nDim_;
    const int irNdim = ir*AlgTraits::nDim_;

    // zero out; scalars that prevail over all components
    double rhoNm1Scs = 0.0;
    double rhoNScs = 0.0;
    double rhoNp1Scs = 0.0;
    double dFdxCont = 0.0;
    double divU = 0.0;

    // zero out vectors that prevail over all components of k
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      v_rhoVrtmScs_(i) = 0.0;
      v_dpdxScs_(i) = 0.0;
    }
    
    // determine scs values of interest
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {

      // save off shape function
      const double r = v_shape_function_(ip,ic);
  
      // time term, density
      rhoNm1Scs += r*v_rhoNm1(ic);
      rhoNScs += r*v_rhoN(ic);
      rhoNp1Scs += r*v_rhoNp1(ic);
      
      // compute scs derivatives and flux derivative
      const double pIC = v_pressure(ic);
      const double rhoIC = v_rhoNp1(ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        const double dnj = v_dndx(ip,ic,j);
        const double vrtmj = v_velocityRTM(ic,j);
        v_rhoVrtmScs_(j) += r*rhoIC*vrtmj;
        divU += r*v_Gju(ic,j,j);
        dFdxCont += rhoIC*vrtmj*dnj;
        v_dpdxScs_(j) += pIC*dnj;
      }
    }
    
    // full continuity residual (constant for all component k)
    const double contRes = (gamma1_*rhoNp1Scs + gamma2_*rhoNScs + gamma3_*rhoNm1Scs)/dt_ + dFdxCont;
    
    // assemble each component
    for ( int k = 0; k < AlgTraits::nDim_; ++k ) {

      const int indexL = ilNdim + k;
      const int indexR = irNdim + k;
      
      const int rowL = indexL*AlgTraits::nodesPerElement_*AlgTraits::nDim_;
      const int rowR = indexR*AlgTraits::nodesPerElement_*AlgTraits::nDim_;

      // zero out residual_k and interpolated velocity_k to scs
      double dFdxkAdv = 0.0;
      double dFdxkDiff = 0.0;
      double ukNm1Scs = 0.0;
      double ukNScs = 0.0;
      double ukNp1Scs = 0.0;

      // zero out vector of local derivatives (duk/dxj)
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        v_dukdxScs_(j) = 0.0;
      }
    
      // determine scs values of interest
      for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
        
        // save off shape function
        const double r = v_shape_function_(ip,ic);
           
        // save off velocityUnp1 for component k
        const double ukNp1 = v_uNp1(ic,k);

        // interpolate all velocity states
        ukNm1Scs += r*v_uNm1(ic,k);
        ukNScs += r*v_uN(ic,k);
        ukNp1Scs += r*ukNp1;
    
        // compute scs derivatives and flux derivative (adv/diff)
        const double rhoIC = v_rhoNp1(ic);
        const double viscIC = v_viscosity(ic);
        for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
          const double dnj = v_dndx(ip,ic,j);
          const double vrtmj = v_velocityRTM(ic,j);
          v_dukdxScs_(j) += ukNp1*dnj;
          dFdxkAdv += rhoIC*vrtmj*ukNp1*dnj;
          dFdxkDiff += viscIC*(v_Gju(ic,k,j) + v_Gju(ic,j,k) - 2.0/3.0*divU*v_kd_(k,j)*includeDivU_)*dnj;
        }
      }
      
      // compute residual for NSO; linearized first
      double residualAlt = dFdxkAdv - ukNp1Scs*dFdxCont;
      for ( int j = 0; j < AlgTraits::nDim_; ++j )
        residualAlt -= v_rhoVrtmScs_(j)*v_dukdxScs_(j);
      
      // compute residual for NSO; pde-based second
      const double time = (gamma1_*rhoNp1Scs*ukNp1Scs + gamma2_*rhoNScs*ukNScs + gamma3_*rhoNm1Scs*ukNm1Scs)/dt_;
      const double residualPde = time + dFdxkAdv - dFdxkDiff + v_dpdxScs_[k] - contRes*ukNp1Scs*nonConservedForm_;

      // final form
      const double residual = residualAlt*altResFac_ + residualPde*om_altResFac_;

      // denominator for nu as well as terms for "upwind" nu
      double gUpperMagGradQ = 0.0;
      double rhoVrtmiGLowerRhoVrtmj = 0.0;
      for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
        const double duidxScs = v_dukdxScs_(i);
        const double rhoVrtmi = v_rhoVrtmScs_(i);
        for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
          gUpperMagGradQ += duidxScs*v_gijUpper(ip,i,j)*v_dukdxScs_(j);
          rhoVrtmiGLowerRhoVrtmj += rhoVrtmi*v_gijLower(ip,i,j)*v_rhoVrtmScs_(j);
        }
      }
      
      // construct nu from residual
      const double nuResidual = std::sqrt((residual*residual)/(gUpperMagGradQ+small_));
      
      // construct nu from first-order-like approach; SNL-internal write-up (eq 209)
      // for now, only include advection as full set of terms is too diffuse
      const double nuFirstOrder = std::sqrt(rhoVrtmiGLowerRhoVrtmj);
      
      // limit based on first order; Cupw_ is a fudge factor similar to Guermond's approach
      const double nu = std::min(Cupw_*nuFirstOrder, nuResidual);
      
      double gijFac = 0.0;
      for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {

        // save off shape function
        const double r = v_shape_function_(ip,ic);

        // find the row
        const int icNdim = ic*AlgTraits::nDim_;

        const int rLkC_k = rowL+icNdim+k;
        const int rRkC_k = rowR+icNdim+k;

        // save of some variables
        const double ukNp1 = v_uNp1(ic,k);
        
        // NSO diffusion-like term; -nu*gUpper*dQ/dxj*ai (residual below)
        double lhsfac = 0.0;
        for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
          const double axi = v_scs_areav(ip,i);
          for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
            const double dnxj = v_dndx(ip,ic,j);
            const double fac = v_gijUpper(ip,i,j)*dnxj*axi;
            const double facGj = r*v_gijUpper(ip,i,j)*v_Gju(ic,k,j)*axi;
            gijFac += fac*ukNp1 - facGj*fourthFac_;
            lhsfac += -fac;
          }
        }
        
        // no coupling between components
        lhs[rLkC_k] += nu*lhsfac;
        lhs[rRkC_k] -= nu*lhsfac;
      }
      
      // residual; left and right
      const double residualNSO = -nu*gijFac;
      rhs[indexL] -= residualNSO;
      rhs[indexR] += residualNSO;
    }
  }
}

INSTANTIATE_SUPPLEMENTAL_ALGORITHM(MomentumNSOElemSuppAlg);
  
} // namespace nalu
} // namespace Sierra
