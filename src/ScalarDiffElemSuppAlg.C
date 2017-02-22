/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ScalarDiffElemSuppAlg.h>
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
#include <stk_mesh/base/BulkData.hpp>
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
// ScalarDiffElemSuppAlg - CVFEM scalar diffusion kernel
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
ScalarDiffElemSuppAlg<AlgTraits>::ScalarDiffElemSuppAlg(
  Realm &realm,
  ScalarFieldType *scalarQ,
  ScalarFieldType *diffFluxCoeff,
  ElemDataRequests& dataPreReqs)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    scalarQ_(scalarQ),
    diffFluxCoeff_(diffFluxCoeff),
    coordinates_(NULL),
    lrscv_(realm.get_surface_master_element(AlgTraits::topo_)->adjacentNodes()),
    ws_shape_function_("ws_shape_function", AlgTraits::numScsIp_, AlgTraits::nodesPerElement_)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // compute shape function; do we want to push this to dataPreReqs?
  MasterElement *meSCS = realm.get_surface_master_element(AlgTraits::topo_);
  meSCS->shape_fcn(&ws_shape_function_(0,0));
  
  // add master elements
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields and data
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
ScalarDiffElemSuppAlg<AlgTraits>::element_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity element,
  ScratchViews& scratchViews)
{
  SharedMemView<double*>& scalarQ_view = scratchViews.get_scratch_view_1D(*scalarQ_);
  SharedMemView<double*>& diffFluxCoeff_view = scratchViews.get_scratch_view_1D(*diffFluxCoeff_);

  SharedMemView<double**>& scs_areav = scratchViews.scs_areav;
  SharedMemView<double***>& dndx = scratchViews.dndx;

  // start the assembly
  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {
    
    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];
    
    // corresponding matrix rows
    const int rowL = il*AlgTraits::nodesPerElement_;
    const int rowR = ir*AlgTraits::nodesPerElement_;

    // compute ip property
    double diffFluxCoeffIp = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const double r = ws_shape_function_(ip,ic);
      diffFluxCoeffIp += r*diffFluxCoeff_view(ic);
    }

    // assemble to rhs and lhs
    double qDiff = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {      
      double lhsfacDiff = 0.0;
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        lhsfacDiff += -diffFluxCoeffIp*dndx(ip,ic,j)*scs_areav(ip,j);
      }
      qDiff += lhsfacDiff*scalarQ_view(ic);
      
      // lhs; il then ir
      lhs[rowL+ic] += lhsfacDiff;
      lhs[rowR+ic] -= lhsfacDiff;
    }
    
    // rhs; il then ir
    rhs[il] -= qDiff;
    rhs[ir] += qDiff;
  }
}

INSTANTIATE_SUPPLEMENTAL_ALGORITHM(ScalarDiffElemSuppAlg);
  
} // namespace nalu
} // namespace Sierra
