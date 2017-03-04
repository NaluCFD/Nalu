/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <MomentumBuoyancySrcElemSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SolutionOptions.h>

// master element
#include <master_element/MasterElement.h>

#include <BuildTemplates.h>
#include <ScratchViews.h>

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
// MomentumBuoyancySrcElemSuppAlg - CMM for momentum buoyancy (u-dof)
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
template<typename AlgTraits>
MomentumBuoyancySrcElemSuppAlg<AlgTraits>::MomentumBuoyancySrcElemSuppAlg(
  Realm &realm,
  ElemDataRequests& dataPreReqs)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    densityNp1_(NULL),
    coordinates_(NULL),
    rhoRef_(0.0),
    useShifted_(false),
    ipNodeMap_(realm.get_volume_master_element(AlgTraits::topo_)->ipNodeMap())
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  
  for (int j=0; j < AlgTraits::nDim_; j++)
    gravity_(j) = realm_.solutionOptions_->gravity_[j];
  rhoRef_ = realm_.solutionOptions_->referenceDensity_;

  MasterElement* meSCV = realm.get_volume_master_element(AlgTraits::topo_);

  meSCV->shape_fcn(&v_shape_function_(0,0));

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_gathered_nodal_field(*coordinates_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_, 1);
  dataPreReqs.add_master_element_call(SCV_VOLUME);
}

//--------------------------------------------------------------------------
//-------- elem_execute ----------------------------------------------------
//--------------------------------------------------------------------------
template<typename AlgTraits>
void
MomentumBuoyancySrcElemSuppAlg<AlgTraits>::element_execute(
  double* /* lhs */,
  double* rhs,
  stk::mesh::Entity /*element*/,
  ScratchViews& scratchViews)
{
  SharedMemView<double*>& v_densityNp1 = scratchViews.get_scratch_view_1D(*densityNp1_);
  SharedMemView<double*>& v_scv_volume = scratchViews.scv_volume;

  for (int ip=0; ip < AlgTraits::numScvIp_; ++ip) {
    const int nearestNode = ipNodeMap_[ip];
    double rhoNp1 = 0.0;

    for (int ic=0; ic < AlgTraits::nodesPerElement_; ++ic) {
      const double r = v_shape_function_(ip, ic);
      rhoNp1 += r * v_densityNp1(ic);
    }

    // Compute RHS
    const double scV = v_scv_volume(ip);
    const int nnNdim = nearestNode * AlgTraits::nDim_;
    const double fac = (rhoNp1 - rhoRef_) * scV;
    for (int j=0; j < AlgTraits::nDim_; j++) {
      rhs[nnNdim + j] += fac * gravity_(j);
    }

    // No LHS contributions
  }
}

INSTANTIATE_SUPPLEMENTAL_ALGORITHM(MomentumBuoyancySrcElemSuppAlg);

} // namespace nalu
} // namespace Sierra
