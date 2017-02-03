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

// manage supplemental algorithms that are templated
#include <BuildTemplates.h>

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
  const stk::topology &theTopo)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    scalarQ_(scalarQ),
    diffFluxCoeff_(diffFluxCoeff),
    coordinates_(NULL),
    meSCS_(realm.get_surface_master_element(theTopo)),
    lrscv_(meSCS_->adjacentNodes()),
    ws_scs_areav_("ws_scs_areav", AlgTraits::numScsIp_, AlgTraits::nDim_),
    ws_dndx_("ws_dndx", AlgTraits::numScsIp_, AlgTraits::nodesPerElement_, AlgTraits::nDim_),
    ws_deriv_("ws_deriv", AlgTraits::nDim_*AlgTraits::numScsIp_*AlgTraits::nodesPerElement_ ),
    ws_det_j_("ws_det_j", AlgTraits::numScsIp_),
    ws_shape_function_("ws_shape_function", AlgTraits::numScsIp_, AlgTraits::nodesPerElement_),
    ws_scalarQ_("ws_scalarQ", AlgTraits::nodesPerElement_),
    ws_diffFluxCoeff_("ws_diffFluxCoeff", AlgTraits::nodesPerElement_),
    ws_coordinates_("ws_coordinates", AlgTraits::nodesPerElement_, AlgTraits::nDim_)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // compute shape function
  meSCS_->shape_fcn(&ws_shape_function_(0,0));
}

//--------------------------------------------------------------------------
//-------- element_execute -------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
void
ScalarDiffElemSuppAlg<AlgTraits>::element_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity element)
{
  // gather
  stk::mesh::Entity const *  node_rels = bulkData_->begin_nodes(element);
  int num_nodes = bulkData_->num_nodes(element);

  // sanity check on num nodes
  ThrowAssert( num_nodes == AlgTraits::nodesPerElement_ );
  for ( int ni = 0; ni < num_nodes; ++ni ) {
    stk::mesh::Entity node = node_rels[ni];
    
    // gather scalars
    ws_scalarQ_(ni) = *stk::mesh::field_data(*scalarQ_, node);
    ws_diffFluxCoeff_(ni) = *stk::mesh::field_data(*diffFluxCoeff_, node );

    // pointers to real data
    const double * coords = stk::mesh::field_data(*coordinates_, node );

    // gather vectors
    for ( int j=0; j < AlgTraits::nDim_; ++j ) {
      ws_coordinates_(ni,j) = coords[j];
    }
  }

  // compute geometry and dndx
  double scs_error = 0.0;
  meSCS_->determinant(1, &ws_coordinates_(0,0), &ws_scs_areav_(0,0), &scs_error);
  meSCS_->grad_op(1, &ws_coordinates_(0,0), &ws_dndx_(0,0,0), &ws_deriv_(0), &ws_det_j_(0), &scs_error);
  
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
      diffFluxCoeffIp += r*ws_diffFluxCoeff_(ic);
    }

    // assemble to rhs and lhs
    double qDiff = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {      
      double lhsfacDiff = 0.0;
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        lhsfacDiff += -diffFluxCoeffIp*ws_dndx_(ip,ic,j)*ws_scs_areav_(ip,j);
      }
      qDiff += lhsfacDiff*ws_scalarQ_(ic);
      
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
