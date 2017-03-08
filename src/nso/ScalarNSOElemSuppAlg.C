/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <nso/ScalarNSOElemSuppAlg.h>
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
// ScalarNSOElemSuppAlg - NSO for scalar equation
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
ScalarNSOElemSuppAlg<AlgTraits>::ScalarNSOElemSuppAlg(
  Realm &realm,
  ScalarFieldType *scalarQ,
  VectorFieldType *Gjq,
  ScalarFieldType *diffFluxCoeff,
  const double fourthFac,
  const double altResFac,
  ElemDataRequests& dataPreReqs)
  : SupplementalAlgorithm(realm),
    scalarQNm1_(NULL),
    scalarQN_(NULL),
    scalarQNp1_(NULL),
    densityNm1_(NULL),
    densityN_(NULL),
    densityNp1_(NULL),
    diffFluxCoeff_(diffFluxCoeff),
    velocityRTM_(NULL),
    Gjq_(Gjq),
    coordinates_(NULL),
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
    nonConservedForm_(0.0)
{
  // save off fields; for non-BDF2 gather in state N for Nm1 (gamma3_ will be zero)
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  scalarQNm1_ = realm_.number_of_states() == 2 ? &(scalarQ->field_of_state(stk::mesh::StateN)) : &(scalarQ->field_of_state(stk::mesh::StateNM1));
  scalarQN_ = &(scalarQ->field_of_state(stk::mesh::StateN));
  scalarQNp1_ = &(scalarQ->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNm1_ = realm_.number_of_states() == 2 ? &(density->field_of_state(stk::mesh::StateN)) : &(density->field_of_state(stk::mesh::StateNM1));
  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
 
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
  dataPreReqs.add_gathered_nodal_field(*coordinates_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityRTM_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*Gjq_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*scalarQNm1_, 1);
  dataPreReqs.add_gathered_nodal_field(*scalarQN_, 1);
  dataPreReqs.add_gathered_nodal_field(*scalarQNp1_, 1);
  
  dataPreReqs.add_gathered_nodal_field(*densityNm1_,1);
  dataPreReqs.add_gathered_nodal_field(*densityN_,1);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_,1);
  dataPreReqs.add_gathered_nodal_field(*diffFluxCoeff_,1);
  
  // master element data
  dataPreReqs.add_master_element_call(SCS_AREAV);
  dataPreReqs.add_master_element_call(SCS_GRAD_OP);
  dataPreReqs.add_master_element_call(SCS_GIJ); 
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
void
ScalarNSOElemSuppAlg<AlgTraits>::setup()
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
ScalarNSOElemSuppAlg<AlgTraits>::element_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity element,
  ScratchViews& scratchViews)
{
  SharedMemView<double**>& v_Gjq = scratchViews.get_scratch_view_2D(*Gjq_);
  SharedMemView<double**>& v_velocityRTM = scratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<double*>& v_qNm1 = scratchViews.get_scratch_view_1D(*scalarQNm1_);
  SharedMemView<double*>& v_qN = scratchViews.get_scratch_view_1D(*scalarQN_);
  SharedMemView<double*>& v_qNp1 = scratchViews.get_scratch_view_1D(*scalarQNp1_);
  SharedMemView<double*>& v_rhoNm1 = scratchViews.get_scratch_view_1D(*densityNm1_);
  SharedMemView<double*>& v_rhoN = scratchViews.get_scratch_view_1D(*densityN_);
  SharedMemView<double*>& v_rhoNp1 = scratchViews.get_scratch_view_1D(*densityNp1_);
  SharedMemView<double*>& v_diffFluxCoeff = scratchViews.get_scratch_view_1D(*diffFluxCoeff_);
 
  SharedMemView<double**>& v_scs_areav = scratchViews.scs_areav;
  SharedMemView<double***>& v_dndx = scratchViews.dndx;
  SharedMemView<double***>& v_gijUpper = scratchViews.gijUpper;
  SharedMemView<double***>& v_gijLower = scratchViews.gijLower;

  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {

    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];

    // corresponding matrix rows
    const int rowL = il*AlgTraits::nodesPerElement_;
    const int rowR = ir*AlgTraits::nodesPerElement_;
   
    // zero out; scalar
    double qNm1Scs = 0.0;
    double qNScs = 0.0;
    double qNp1Scs = 0.0;
    double rhoNm1Scs = 0.0;
    double rhoNScs = 0.0;
    double rhoNp1Scs = 0.0;
    double dFdxAdv = 0.0;
    double dFdxDiff = 0.0;
    double dFdxCont = 0.0;
   
    // zero out vector
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      v_dqdxScs_(i) = 0.0;
      v_rhoVrtmScs_(i) = 0.0;
    }
    
    // determine scs values of interest
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      // save off shape function
      const double r = v_shape_function_(ip,ic);

      // time term; scalar q
      qNm1Scs += r*v_qNm1(ic);
      qNScs += r*v_qN(ic);
      qNp1Scs += r*v_qNp1(ic);

      // time term, density
      rhoNm1Scs += r*v_rhoNm1(ic);
      rhoNScs += r*v_rhoN(ic);
      rhoNp1Scs += r*v_rhoNp1(ic);

      // compute scs derivatives and flux derivative
      const double qIC = v_qNp1(ic);
      const double rhoIC = v_rhoNp1(ic);
      const double diffFluxCoeffIC = v_diffFluxCoeff(ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        const double dnj = v_dndx(ip,ic,j);
        const double vrtmj = v_velocityRTM(ic,j);
        v_dqdxScs_(j) += qIC*dnj;
        v_rhoVrtmScs_(j) += r*rhoIC*vrtmj;
        dFdxAdv += rhoIC*vrtmj*qIC*dnj;
        dFdxDiff += diffFluxCoeffIC*v_Gjq(ic,j)*dnj;
        dFdxCont += rhoIC*vrtmj*dnj;
      }
    }
    
    // full continuity residual
    const double contRes = (gamma1_*rhoNp1Scs + gamma2_*rhoNScs + gamma3_*rhoNm1Scs)/dt_ + dFdxCont;

    // compute residual for NSO; linearized first
    double residualAlt = dFdxAdv - qNp1Scs*dFdxCont;
    for ( int j = 0; j < AlgTraits::nDim_; ++j )
      residualAlt -= v_rhoVrtmScs_(j)*v_dqdxScs_(j);
    
    // compute residual for NSO; pde-based second
    const double time = (gamma1_*rhoNp1Scs*qNp1Scs + gamma2_*rhoNScs*qNScs + gamma3_*rhoNm1Scs*qNm1Scs)/dt_;
    const double residualPde = time + dFdxAdv - dFdxDiff - contRes*qNp1Scs*nonConservedForm_;

    // final form
    const double residual = residualAlt*altResFac_ + residualPde*om_altResFac_;

    // denominator for nu as well as terms for "upwind" nu
    double gUpperMagGradQ = 0.0;
    double rhoVrtmiGLowerRhoVrtmj = 0.0;
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      const double dqdxScsi = v_dqdxScs_(i);
      const double rhoVrtmi = v_rhoVrtmScs_(i);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        gUpperMagGradQ += dqdxScsi*v_gijUpper(ip,i,j)*v_dqdxScs_[j];
        rhoVrtmiGLowerRhoVrtmj += rhoVrtmi*v_gijLower(ip,i,j)*v_rhoVrtmScs_[j];
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
      
      // save of some variables
      const double qIC = v_qNp1(ic);
      
      // NSO diffusion-like term; -nu*gUpper*(dQ/dxj - Gjq)*ai (residual below)
      double lhsfac = 0.0;
      for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
        const double axi = v_scs_areav(ip,i);
        for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
          const double dnxj = v_dndx(ip,ic,j);
          const double fac = v_gijUpper(ip,i,j)*dnxj*axi;
          const double facGj = r*v_gijUpper(ip,i,j)*v_Gjq(ic,j)*axi;
          gijFac += fac*qIC - facGj*fourthFac_;
          lhsfac += -fac;
        }
      }
      
      lhs[rowL+ic] += nu*lhsfac;
      lhs[rowR+ic] -= nu*lhsfac;
    }
    
    // residual; left and right
    const double residualNSO = -nu*gijFac;
    rhs[il] -= residualNSO;
    rhs[ir] += residualNSO;
  }      
}

INSTANTIATE_SUPPLEMENTAL_ALGORITHM(ScalarNSOElemSuppAlg);
  
} // namespace nalu
} // namespace Sierra
