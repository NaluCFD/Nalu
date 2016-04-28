/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <MomentumActuatorLineSrcNodeSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>

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
// MomentumActuatorLineSrcNodeSuppAlg - actuator line drag source term
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MomentumActuatorLineSrcNodeSuppAlg::MomentumActuatorLineSrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    actuatorLineSrc_(NULL),
    actuatorLineSrcLHS_(NULL),
    dualNodalVolume_(NULL),
    nDim_(1)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  actuatorLineSrc_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "actuator_line_source");
  actuatorLineSrcLHS_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "actuator_line_source_lhs");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  nDim_ = meta_data.spatial_dimension();
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumActuatorLineSrcNodeSuppAlg::setup()
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumActuatorLineSrcNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // single point quadrature
  const double *src = stk::mesh::field_data(*actuatorLineSrc_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node);
  const double srcLHS = *stk::mesh::field_data(*actuatorLineSrcLHS_, node);
 
  const int nDim = nDim_;
  for ( int i = 0; i < nDim; ++i ) {
    rhs[i] += src[i]*dualVolume;
    const int row = i*nDim;
    lhs[row+i] += srcLHS;
  }
}

} // namespace nalu
} // namespace Sierra
