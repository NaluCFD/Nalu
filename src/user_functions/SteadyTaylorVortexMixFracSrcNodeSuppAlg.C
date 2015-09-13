/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/SteadyTaylorVortexMixFracSrcNodeSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <TimeIntegrator.h>

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
// SteadyTaylorVortexMixFracSrcNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SteadyTaylorVortexMixFracSrcNodeSuppAlg::SteadyTaylorVortexMixFracSrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    coordinates_(NULL),
    dualNodalVolume_(NULL),
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
    pi_(acos(-1.0))
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
SteadyTaylorVortexMixFracSrcNodeSuppAlg::setup()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
SteadyTaylorVortexMixFracSrcNodeSuppAlg::node_execute(
  double */*lhs*/,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix
  const double *coords = stk::mesh::field_data(*coordinates_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double x = coords[0];
  const double y = coords[1];

  const double src = amf_ * pi_ / Sc_ * (cos(a_ * pi_ * x) * sin(a_ * pi_ * y) * sin(amf_ * pi_ * x) * sin(amf_ * pi_ * y) * Sc_ + sin(a_ * pi_ * x) * cos(a_ * pi_ * y) * cos(amf_ * pi_ * x) * cos(amf_ * pi_ * y) * Sc_ + 0.2e1 * visc_ * cos(amf_ * pi_ * x) * amf_ * pi_ * sin(amf_ * pi_ * y));

  rhs[0] += src*dualVolume;
}

} // namespace nalu
} // namespace Sierra
