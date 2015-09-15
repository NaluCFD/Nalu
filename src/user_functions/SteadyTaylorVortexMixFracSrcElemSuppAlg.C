/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/SteadyTaylorVortexMixFracSrcElemSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// SteadyTaylorVortexMixFracSrcElemSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SteadyTaylorVortexMixFracSrcElemSuppAlg::SteadyTaylorVortexMixFracSrcElemSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    coordinates_(NULL),
    nDim_(realm_.spatialDimension_),
    rhoP_(1.0),
    rhoS_(1.0),
    unot_(1.0),
    vnot_(1.0),
    znot_(1.0),
    pnot_(1.0),
    visc_(0.001),
    a_(20.0),
    amf_(10.0),
    Sc_(0.9),
    pi_(acos(-1.0)),
    useShifted_(false)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
 
  // scratch vecs
  scvCoords_.resize(nDim_);
}

//--------------------------------------------------------------------------
//-------- elem_resize -----------------------------------------------------
//--------------------------------------------------------------------------
void
SteadyTaylorVortexMixFracSrcElemSuppAlg::elem_resize(
  MasterElement */*meSCS*/,
  MasterElement *meSCV)
{
  const int nodesPerElement = meSCV->nodesPerElement_;
  const int numScvIp = meSCV->numIntPoints_;
  ws_shape_function_.resize(numScvIp*nodesPerElement);
  ws_coordinates_.resize(nDim_*nodesPerElement);
  ws_scv_volume_.resize(numScvIp);

  // compute shape function
  const bool doIt = true;
  if ( doIt ) {
  if ( useShifted_ )
    meSCV->shifted_shape_fcn(&ws_shape_function_[0]);
  else
    meSCV->shape_fcn(&ws_shape_function_[0]);
  }
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
SteadyTaylorVortexMixFracSrcElemSuppAlg::setup()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- elem_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
SteadyTaylorVortexMixFracSrcElemSuppAlg::elem_execute(
  double */*lhs*/,
  double *rhs,
  stk::mesh::Entity element,
  MasterElement */*meSCS*/,
  MasterElement *meSCV)
{
  // pointer to ME methods
  const int *ipNodeMap = meSCV->ipNodeMap();
  const int nodesPerElement = meSCV->nodesPerElement_;
  const int numScvIp = meSCV->numIntPoints_;
  
  // gather
  stk::mesh::Entity const *  node_rels = bulkData_->begin_nodes(element);
  int num_nodes = bulkData_->num_nodes(element);
  
  // sanity check on num nodes
  ThrowAssert( num_nodes == nodesPerElement );
  
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
  meSCV->determinant(1, &ws_coordinates_[0], &ws_scv_volume_[0], &scv_error);

  for ( int ip = 0; ip < numScvIp; ++ip ) {
    
    // nearest node to ip
    const int nearestNode = ipNodeMap[ip];
    
    // zero out
    for ( int j =0; j < nDim_; ++j )
      scvCoords_[j] = 0.0;
       
    const int offSet = ip*nodesPerElement;
    for ( int ic = 0; ic < nodesPerElement; ++ic ) {
      const double r = ws_shape_function_[offSet+ic];
      for ( int j = 0; j < nDim_; ++j )
        scvCoords_[j] += r*ws_coordinates_[ic*nDim_+j];
    }
    const double x = scvCoords_[0];
    const double y = scvCoords_[1];
    
    const double src = amf_ * pi_ / Sc_ * (cos(a_ * pi_ * x) * sin(a_ * pi_ * y) * sin(amf_ * pi_ * x) * sin(amf_ * pi_ * y) * Sc_ + sin(a_ * pi_ * x) * cos(a_ * pi_ * y) * cos(amf_ * pi_ * x) * cos(amf_ * pi_ * y) * Sc_ + 0.2e1 * visc_ * cos(amf_ * pi_ * x) * amf_ * pi_ * sin(amf_ * pi_ * y));

    rhs[nearestNode] += src*ws_coordinates_[ip];      
  }
}

} // namespace nalu
} // namespace Sierra
