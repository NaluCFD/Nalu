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
// SteadyThermal3dContactSrcElemSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
SteadyThermal3dContactSrcElemSuppAlg<AlgTraits>::SteadyThermal3dContactSrcElemSuppAlg(
  Realm &realm,
  ElemDataRequests& dataPreReqs)
  : SupplementalAlgorithm(realm),
    coordinates_(NULL),
    ipNodeMap_(realm.get_volume_master_element(AlgTraits::topo_)->ipNodeMap()),
    a_(1.0),
    k_(1.0),
    pi_(std::acos(-1.0))
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
 
  // compute shape function; possibly push this to dataPreReqs?
  MasterElement *meSCV = realm.get_volume_master_element(AlgTraits::topo_);
  meSCV->shape_fcn(&v_shape_function_(0,0));

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_gathered_nodal_field(*coordinates_, AlgTraits::nDim_);
  dataPreReqs.add_master_element_call(SCV_VOLUME);
}

//--------------------------------------------------------------------------
//-------- element_execute -------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
void
SteadyThermal3dContactSrcElemSuppAlg<AlgTraits>::element_execute(
  double */*lhs*/,
  double *rhs,
  stk::mesh::Entity element,
  ScratchViews& scratchViews)
{
  SharedMemView<double**>& v_coordinates = scratchViews.get_scratch_view_2D(*coordinates_);
  SharedMemView<double*>& v_scv_volume = scratchViews.scv_volume;

  // interpolate to ips and evaluate source
  for ( int ip = 0; ip < AlgTraits::numScvIp_; ++ip ) {
    
    // nearest node to ip
    const int nearestNode = ipNodeMap_[ip];
    
    // zero out
    for ( int j =0; j < AlgTraits::nDim_; ++j )
      v_scvCoords_(j) = 0.0;
    
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const double r = v_shape_function_(ip,ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j )
        v_scvCoords_(j) += r*v_coordinates(ic,j);
    }
    const double x = v_scvCoords_(0);
    const double y = v_scvCoords_(1);
    const double z = v_scvCoords_(2);
    rhs[nearestNode] += k_/4.0*(2.0*a_*pi_)*(2.0*a_*pi_)*(cos(2.0*a_*pi_*x) 
                                                          + cos(2.0*a_*pi_*y) 
                                                          + cos(2.0*a_*pi_*z))*v_scv_volume(ip);
  }
}

INSTANTIATE_SUPPLEMENTAL_ALGORITHM(SteadyThermal3dContactSrcElemSuppAlg);

} // namespace nalu
} // namespace Sierra
