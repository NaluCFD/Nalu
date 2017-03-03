/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ScalarAdvDiffElemSuppAlg.h>
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
// ScalarAdvDiffElemSuppAlg - CVFEM scalar adv/diffusion kernel
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
ScalarAdvDiffElemSuppAlg<AlgTraits>::ScalarAdvDiffElemSuppAlg(
  Realm &realm,
  ScalarFieldType *scalarQ,
  ScalarFieldType *diffFluxCoeff,
  ElemDataRequests& dataPreReqs)
  : SupplementalAlgorithm(realm),
    scalarQ_(scalarQ),
    diffFluxCoeff_(diffFluxCoeff),
    coordinates_(NULL),
    massFlowRate_(NULL),
    lrscv_(realm.get_surface_master_element(AlgTraits::topo_)->adjacentNodes())
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  massFlowRate_ = meta_data.get_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "mass_flow_rate_scs");

  // compute shape function; do we want to push this to dataPreReqs?
  MasterElement *meSCS = realm.get_surface_master_element(AlgTraits::topo_);
  meSCS->shape_fcn(&v_shape_function_(0,0));
  
  // add master elements
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields and data; mdot not gathered
  dataPreReqs.add_gathered_nodal_field(*coordinates_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*scalarQ, 1);
  dataPreReqs.add_gathered_nodal_field(*diffFluxCoeff, 1);
  dataPreReqs.add_master_element_call(SCS_AREAV);
  dataPreReqs.add_master_element_call(SCS_GRAD_OP);
}

//--------------------------------------------------------------------------
//-------- element_execute -------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
void
ScalarAdvDiffElemSuppAlg<AlgTraits>::element_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity element,
  ScratchViews& scratchViews)
{
  SharedMemView<double*>& v_scalarQ = scratchViews.get_scratch_view_1D(*scalarQ_);
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
    
    // corresponding matrix rows
    const int rowL = il*AlgTraits::nodesPerElement_;
    const int rowR = ir*AlgTraits::nodesPerElement_;

    // save off mdot
    const double tmdot = mdot[ip];

    // compute ip property and 
    double diffFluxCoeffIp = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const double r = v_shape_function_(ip,ic);
      diffFluxCoeffIp += r*v_diffFluxCoeff(ic);
    }
      
    // advection and diffusion 
    double qAdv = 0.0;
    double qDiff = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {      

      // advection
      const double r = v_shape_function_(ip,ic);
      const double lhsfacAdv = r*tmdot;
      qAdv += lhsfacAdv*v_scalarQ(ic);

      // diffusion
      double lhsfacDiff = 0.0;
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        lhsfacDiff += -diffFluxCoeffIp*v_dndx(ip,ic,j)*v_scs_areav(ip,j);
      }
      qDiff += lhsfacDiff*v_scalarQ(ic);
      
      // lhs; il then ir
      lhs[rowL+ic] += lhsfacAdv + lhsfacDiff;
      lhs[rowR+ic] -= lhsfacAdv + lhsfacDiff;
    }
    
    // rhs; il then ir
    rhs[il] -= qAdv + qDiff;
    rhs[ir] += qAdv + qDiff;
  }
}

INSTANTIATE_SUPPLEMENTAL_ALGORITHM(ScalarAdvDiffElemSuppAlg);
  
} // namespace nalu
} // namespace Sierra
