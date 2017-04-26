/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ScalarUpwAdvDiffElemSuppAlg.h>
#include <EquationSystem.h>
#include <SupplementalAlgorithm.h>

#include <FieldTypeDef.h>
#include <Realm.h>
#include <PecletFunction.h>
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
// ScalarUpwAdvDiffElemSuppAlg - CVFEM scalar adv/diffusion kernel
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
ScalarUpwAdvDiffElemSuppAlg<AlgTraits>::ScalarUpwAdvDiffElemSuppAlg(
  Realm &realm, 
  EquationSystem *eqSystem,
  ScalarFieldType *scalarQ,
  VectorFieldType *Gjq,
  ScalarFieldType *diffFluxCoeff,
  ElemDataRequests& dataPreReqs)
  : SupplementalAlgorithm(realm),
    scalarQ_(scalarQ),
    Gjq_(Gjq),
    diffFluxCoeff_(diffFluxCoeff),
    velocityRTM_(NULL),
    coordinates_(NULL),
    density_(NULL),
    massFlowRate_(NULL),
    lrscv_(sierra::nalu::get_surface_master_element(AlgTraits::topo_)->adjacentNodes()),
    dofName_(scalarQ_->name()),
    alpha_(realm.get_alpha_factor(dofName_)),
    alphaUpw_(realm.get_alpha_upw_factor(dofName_)),
    hoUpwind_(realm.get_upw_factor(dofName_)),
    useLimiter_(realm.primitive_uses_limiter(dofName_)),
    om_alpha_(1.0-alpha_),
    om_alphaUpw_(1.0-alphaUpw_),
    small_(1.0e-16),
    pecletFunction_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  if ( realm_.does_mesh_move() )
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");
  else 
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  massFlowRate_ = meta_data.get_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "mass_flow_rate_scs");

  // create the peclet blending function
  pecletFunction_ = eqSystem->create_peclet_function(scalarQ_->name());

  // compute shape function; do we want to push this to dataPreReqs?
  MasterElement *meSCS = sierra::nalu::get_surface_master_element(AlgTraits::topo_);
  meSCS->shape_fcn(&v_shape_function_(0,0));
  
  // add master elements
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields and data; mdot not gathered
  dataPreReqs.add_gathered_nodal_field(*velocityRTM_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*coordinates_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*Gjq, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*scalarQ, 1);
  dataPreReqs.add_gathered_nodal_field(*density_, 1);
  dataPreReqs.add_gathered_nodal_field(*diffFluxCoeff, 1);
  dataPreReqs.add_master_element_call(SCS_AREAV);
  dataPreReqs.add_master_element_call(SCS_GRAD_OP);
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
ScalarUpwAdvDiffElemSuppAlg<AlgTraits>::~ScalarUpwAdvDiffElemSuppAlg()
{
  delete pecletFunction_;
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
void
ScalarUpwAdvDiffElemSuppAlg<AlgTraits>::setup()
{
  alpha_ = realm_.get_alpha_factor(dofName_);
  alphaUpw_ = realm_.get_alpha_upw_factor(dofName_);
  hoUpwind_ = realm_.get_upw_factor(dofName_);
  useLimiter_ = realm_.primitive_uses_limiter(dofName_);

  // one minus flavor..
  om_alpha_ = 1.0-alpha_;
  om_alphaUpw_ = 1.0-alphaUpw_;
}

//--------------------------------------------------------------------------
//-------- element_execute -------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
void
ScalarUpwAdvDiffElemSuppAlg<AlgTraits>::element_execute(
  SharedMemView<double **>& lhs,
  SharedMemView<double *>& rhs,
  stk::mesh::Entity element,
  ScratchViews& scratchViews)
{
  SharedMemView<double**>& v_velocityRTM = scratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<double**>& v_coordinates = scratchViews.get_scratch_view_2D(*coordinates_);
  SharedMemView<double**>& v_Gjq = scratchViews.get_scratch_view_2D(*Gjq_);
  SharedMemView<double*>& v_scalarQ = scratchViews.get_scratch_view_1D(*scalarQ_);
  SharedMemView<double*>& v_density = scratchViews.get_scratch_view_1D(*density_);
  SharedMemView<double*>& v_diffFluxCoeff = scratchViews.get_scratch_view_1D(*diffFluxCoeff_);

  SharedMemView<double**>& v_scs_areav = scratchViews.scs_areav;
  SharedMemView<double***>& v_dndx = scratchViews.dndx;

  // ip data for this element
  const double *mdot = stk::mesh::field_data(*massFlowRate_, element);

  // start the assembly
  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {
    
    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];
    
    // save off mdot
    const double tmdot = mdot[ip];

    // zero out values of interest for this ip
    for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
      v_coordIp_(j) = 0.0;
    }

    // compute ip property and 
    double qIp = 0.0;
    double rhoIp = 0.0;
    double diffFluxCoeffIp = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const double r = v_shape_function_(ip,ic);
      qIp += r*v_scalarQ(ic);
      rhoIp += r*v_density(ic);
      diffFluxCoeffIp += r*v_diffFluxCoeff(ic);
      for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
        v_coordIp_(i) += r*v_coordinates(ic,i);
      }
    }
      
    // Peclet factor; along the edge
    const double diffIp = 0.5*(v_diffFluxCoeff(il)/v_density(il)
                               + v_diffFluxCoeff(ir)/v_density(ir));
    double udotx = 0.0;
    for(int j = 0; j < AlgTraits::nDim_; ++j ) {
      const double dxj = v_coordinates(ir,j) - v_coordinates(il,j);;
      const double uj = 0.5*(v_velocityRTM(il,j) + v_velocityRTM(ir,j));
      udotx += uj*dxj;
    }
    const double pecfac = pecletFunction_->execute(std::abs(udotx)/(diffIp+small_));
    const double om_pecfac = 1.0-pecfac;
    
    // left and right extrapolation
    double dqL = 0.0;
    double dqR = 0.0;
    for(int j = 0; j < AlgTraits::nDim_; ++j ) {
      const double dxjL = v_coordIp_(j) - v_coordinates(il,j);
      const double dxjR = v_coordinates(ir,j) - v_coordIp_(j);
      dqL += dxjL*v_Gjq(il,j);
      dqR += dxjR*v_Gjq(ir,j);
    }
    
    // add limiter if appropriate
    double limitL = 1.0;
    double limitR = 1.0;
    if ( useLimiter_ ) {
      const double dq = v_scalarQ(ir) - v_scalarQ(il);
      const double dqMl = 2.0*2.0*dqL - dq;
      const double dqMr = 2.0*2.0*dqR - dq;
      limitL = van_leer(dqMl, dq);
      limitR = van_leer(dqMr, dq);
    }
    
    // extrapolated; for now limit (along edge is fine)
    const double qIpL = v_scalarQ(il) + dqL*hoUpwind_*limitL;
    const double qIpR = v_scalarQ(ir) - dqR*hoUpwind_*limitR;

    // upwind
    const double qUpwind = (tmdot > 0) ? alphaUpw_*qIpL + om_alphaUpw_*qIp
      : alphaUpw_*qIpR + om_alphaUpw_*qIp;

    // generalized central (2nd and 4th order)
    const double qHatL = alpha_*qIpL + om_alpha_*qIp;
    const double qHatR = alpha_*qIpR + om_alpha_*qIp;
    const double qCds = 0.5*(qHatL + qHatR);
    
    // total advection
    const double aflux = tmdot*(pecfac*qUpwind + om_pecfac*qCds);

    // right hand side; L and R
    rhs(il) -= aflux;
    rhs(ir) += aflux; 

    // upwind advection (includes 4th); left node
    const double alhsfacL = 0.5*(tmdot+std::abs(tmdot))*pecfac*alphaUpw_
      + 0.5*alpha_*om_pecfac*tmdot;
    lhs(il,il) += alhsfacL;
    lhs(ir,il) -= alhsfacL;
    
    // upwind advection; right node
    const double alhsfacR = 0.5*(tmdot-std::abs(tmdot))*pecfac*alphaUpw_
      + 0.5*alpha_*om_pecfac*tmdot;
    lhs(ir,ir) -= alhsfacR;
    lhs(il,ir) += alhsfacR;

    // advection and diffusion 
    double qDiff = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {      

      // shape function
      const double r = v_shape_function_(ip,ic);

      // upwind (il/ir) handled above; collect terms on alpha and alphaUpw
      const double lhsfacAdv = r*tmdot*(pecfac*om_alphaUpw_ + om_pecfac*om_alpha_);
      
      // advection operator lhs; rhs handled above
      lhs(il,ic) += lhsfacAdv;
      lhs(ir,ic) -= lhsfacAdv;

      // diffusion
      double lhsfacDiff = 0.0;
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        lhsfacDiff += -diffFluxCoeffIp*v_dndx(ip,ic,j)*v_scs_areav(ip,j);
      }
      qDiff += lhsfacDiff*v_scalarQ(ic);
      
      // lhs; il then ir
      lhs(il,ic) +=  lhsfacDiff;
      lhs(ir,ic) -= lhsfacDiff;
    }
    
    // rhs; il then ir
    rhs(il) -= qDiff;
    rhs(ir) += qDiff;
  }
}

//--------------------------------------------------------------------------
//-------- van_leer --------------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
double
ScalarUpwAdvDiffElemSuppAlg<AlgTraits>::van_leer(
  const double &dqm,
  const double &dqp)
{
  double limit = (2.0*(dqm*dqp+std::abs(dqm*dqp))) /
    ((dqm+dqp)*(dqm+dqp)+small_);
  return limit;
}

INSTANTIATE_SUPPLEMENTAL_ALGORITHM(ScalarUpwAdvDiffElemSuppAlg);
  
} // namespace nalu
} // namespace Sierra
