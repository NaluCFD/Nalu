/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ScalarGclNodeSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <TimeIntegrator.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
// stk_io
#include <stk_io/StkMeshIoBroker.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// ScalarGclNodeSuppAlg - GCL
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ScalarGclNodeSuppAlg::ScalarGclNodeSuppAlg(
  ScalarFieldType *scalarQNp1,
  Realm &realm)
  : SupplementalAlgorithm(realm),
    scalarQNp1_(NULL),
    densityNp1_(NULL),
    dvdx_(NULL),
    dualNodalVolume_(NULL),
    nDim_(1)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  dvdx_ = meta_data.get_field<GenericFieldType>(stk::topology::NODE_RANK, "dvdx");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  nDim_ = meta_data.spatial_dimension();
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
ScalarGclNodeSuppAlg::setup()
{
  // all set up in constructor
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
ScalarGclNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // rhs -= rho*scalarQ*div(v)*dV
  const double scalarQNp1 = *stk::mesh::field_data(*scalarQNp1_, node );
  const double rhoNp1 = *stk::mesh::field_data(*densityNp1_, node );
  const double *dvdx = stk::mesh::field_data(*dvdx_, node );
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );

  // form d/dxj(vj)
  const int nDim = nDim_;
  double divV = 0.0;
  for ( int j = 0; j < nDim; ++j )
    divV += dvdx[j*nDim+j];
  
  rhs[0] -= rhoNp1*scalarQNp1*divV*dualVolume;
  lhs[0] += 0.0;
}

} // namespace nalu
} // namespace Sierra
