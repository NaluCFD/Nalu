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

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

// topology
#include <stk_topology/topology.hpp>

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
SteadyThermal3dContactSrcElemSuppAlg::SteadyThermal3dContactSrcElemSuppAlg(
  Realm &realm,
  const stk::topology &theTopo)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    coordinates_(NULL),
    meSCV_(realm.get_volume_master_element(theTopo)),
    ipNodeMap_(meSCV_->ipNodeMap()),
    nodesPerElement_(meSCV_->nodesPerElement_),
    numScvIp_(meSCV_->numIntPoints_),
    a_(1.0),
    k_(1.0),
    pi_(std::acos(-1.0)),
    useShifted_(false),
    nDim_(realm_.spatialDimension_),
    evalAtIps_(true)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
 
  // scratch vecs; fixed
  scvCoords_.resize(nDim_);

  ws_shape_function_.resize(numScvIp_*nodesPerElement_);
  ws_coordinates_.resize(nDim_*nodesPerElement_);
  ws_scv_volume_.resize(numScvIp_);
  ws_nodalSrc_.resize(nodesPerElement_);

  // compute shape function
  if ( useShifted_ )
    meSCV_->shifted_shape_fcn(&ws_shape_function_[0]);
  else
    meSCV_->shape_fcn(&ws_shape_function_[0]);
}

//--------------------------------------------------------------------------
//-------- element_execute -------------------------------------------------
//--------------------------------------------------------------------------
void
SteadyThermal3dContactSrcElemSuppAlg::element_execute(
  double */*lhs*/,
  double *rhs,
  stk::mesh::Entity element)
{
  // gather
  stk::mesh::Entity const *  node_rels = bulkData_->begin_nodes(element);
  int num_nodes = bulkData_->num_nodes(element);

  // sanity check on num nodes
  ThrowAssert( num_nodes == nodesPerElement_ );

  for ( int ni = 0; ni < num_nodes; ++ni ) {
    stk::mesh::Entity node = node_rels[ni];
    // pointers to real data
    const double * coords = stk::mesh::field_data(*coordinates_, node );
    // gather vectors
    const int niNdim = ni*nDim_;
    for ( int j=0; j < nDim_; ++j ) {
      ws_coordinates_[niNdim+j] = coords[j];
    }
  }

  // compute geometry
  double scv_error = 0.0;
  meSCV_->determinant(1, &ws_coordinates_[0], &ws_scv_volume_[0], &scv_error);

  // choose a form...
  if ( evalAtIps_ ) {
    // interpolate to ips and evaluate source
    for ( int ip = 0; ip < numScvIp_; ++ip ) {
      
      // nearest node to ip
      const int nearestNode = ipNodeMap_[ip];
      
      // zero out
      for ( int j =0; j < nDim_; ++j )
        scvCoords_[j] = 0.0;
      
      const int offSet = ip*nodesPerElement_;
      for ( int ic = 0; ic < nodesPerElement_; ++ic ) {
        const double r = ws_shape_function_[offSet+ic];
        for ( int j = 0; j < nDim_; ++j )
          scvCoords_[j] += r*ws_coordinates_[ic*nDim_+j];
      }
      const double x = scvCoords_[0];
      const double y = scvCoords_[1];
      const double z = scvCoords_[2];
      rhs[nearestNode] += k_/4.0*(2.0*a_*pi_)*(2.0*a_*pi_)*(cos(2.0*a_*pi_*x) + cos(2.0*a_*pi_*y) + cos(2.0*a_*pi_*z))*ws_scv_volume_[ip];
    }
  }
  else {
    // evaluate source at nodal location
    for ( int ni = 0; ni < nodesPerElement_; ++ni ) {
      const double x = ws_coordinates_[ni*nDim_+0];
      const double y = ws_coordinates_[ni*nDim_+1];
      const double z = ws_coordinates_[ni*nDim_+2];
      ws_nodalSrc_[ni] = k_/4.0*(2.0*a_*pi_)*(2.0*a_*pi_)*(cos(2.0*a_*pi_*x) + cos(2.0*a_*pi_*y) + cos(2.0*a_*pi_*z));
    }

    // interpolate nodal source term to ips and assemble to nodes
    for ( int ip = 0; ip < numScvIp_; ++ip ) {
      
      // nearest node to ip
      const int nearestNode = ipNodeMap_[ip];
      
      double ipSource = 0.0;

      const int offSet = ip*nodesPerElement_;
      for ( int ic = 0; ic < nodesPerElement_; ++ic ) {
        const double r = ws_shape_function_[offSet+ic];
        ipSource += r*ws_nodalSrc_[ic];
      }
      rhs[nearestNode] += ipSource*ws_scv_volume_[ip];
    } 
  }
}

} // namespace nalu
} // namespace Sierra
