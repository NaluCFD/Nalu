/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/SteadyThermal3dContactSrcElemSuppAlg.h>
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
// SteadyThermal3dContactSrcElemSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
SteadyThermal3dContactSrcElemSuppAlg<AlgTraits>::SteadyThermal3dContactSrcElemSuppAlg(
  Realm &realm,
  const stk::topology &theTopo)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    coordinates_(NULL),
    meSCV_(realm.get_volume_master_element(theTopo)),
    ipNodeMap_(meSCV_->ipNodeMap()),
    a_(1.0),
    k_(1.0),
    pi_(std::acos(-1.0)),
    ws_shape_function_("ws_shape_function", AlgTraits::numScvIp_, AlgTraits::nodesPerElement_),
    ws_coordinates_("ws_coordinates", AlgTraits::nodesPerElement_, AlgTraits::nDim_),
    ws_scv_volume_("ws_scv_volume", AlgTraits::numScvIp_),
    ws_scvCoords_("ws_scvCoords", AlgTraits::nDim_)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
 
  // compute shape function
  meSCV_->shape_fcn(&ws_shape_function_(0,0));
}

//--------------------------------------------------------------------------
//-------- element_execute -------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
void
SteadyThermal3dContactSrcElemSuppAlg<AlgTraits>::element_execute(
  double */*lhs*/,
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
    // pointers to real data
    const double * coords = stk::mesh::field_data(*coordinates_, node );
    // gather vectors
    for ( int j=0; j < AlgTraits::nDim_; ++j ) {
      ws_coordinates_(ni,j) = coords[j];
    }
  }

  // compute geometry
  double scv_error = 0.0;
  meSCV_->determinant(1, &ws_coordinates_(0,0), &ws_scv_volume_(0), &scv_error);

  // interpolate to ips and evaluate source
  for ( int ip = 0; ip < AlgTraits::numScvIp_; ++ip ) {
    
    // nearest node to ip
    const int nearestNode = ipNodeMap_[ip];
    
    // zero out
    for ( int j =0; j < AlgTraits::nDim_; ++j )
      ws_scvCoords_(j) = 0.0;
    
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const double r = ws_shape_function_(ip,ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j )
        ws_scvCoords_(j) += r*ws_coordinates_(ic,j);
    }
    const double x = ws_scvCoords_(0);
    const double y = ws_scvCoords_(1);
    const double z = ws_scvCoords_(2);
    rhs[nearestNode] += k_/4.0*(2.0*a_*pi_)*(2.0*a_*pi_)*(cos(2.0*a_*pi_*x) 
                                                          + cos(2.0*a_*pi_*y) 
                                                          + cos(2.0*a_*pi_*z))*ws_scv_volume_[ip];
  }
}

INSTANTIATE_SUPPLEMENTAL_ALGORITHM(SteadyThermal3dContactSrcElemSuppAlg);

} // namespace nalu
} // namespace Sierra
