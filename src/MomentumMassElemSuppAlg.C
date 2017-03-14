/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <MomentumMassElemSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
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
// MomentumMassElemSuppAlg - CMM (BDF2/BE) for momentum equation (u-dof)
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
template<typename AlgTraits>
MomentumMassElemSuppAlg<AlgTraits>::MomentumMassElemSuppAlg(
  Realm &realm,
  ElemDataRequests& dataPreReqs,
   const bool lumpedMass)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    velocityNm1_(NULL),
    velocityN_(NULL),
    velocityNp1_(NULL),
    densityNm1_(NULL),
    densityN_(NULL),
    densityNp1_(NULL),
    Gjp_(NULL),
    coordinates_(NULL),
    dt_(0.0),
    gamma1_(0.0),
    gamma2_(0.0),
    gamma3_(0.0),
    lumpedMass_(lumpedMass),
    ipNodeMap_(realm.get_volume_master_element(AlgTraits::topo_)->ipNodeMap())
{
  // save off fields; shove state N into Nm1 if this is BE
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  VectorFieldType *velocity = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  velocityNm1_ = realm_.number_of_states() == 2 ? &(velocity->field_of_state(stk::mesh::StateN)) : &(velocity->field_of_state(stk::mesh::StateNM1));
  velocityN_ = &(velocity->field_of_state(stk::mesh::StateN));
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNm1_ = realm_.number_of_states() == 2 ? &(density->field_of_state(stk::mesh::StateN)) : &(density->field_of_state(stk::mesh::StateNM1));
  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  Gjp_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  MasterElement* meSCV = realm.get_volume_master_element(AlgTraits::topo_);

  // compute shape function
  if ( lumpedMass_ )
    meSCV->shifted_shape_fcn(&v_shape_function_(0,0));
  else
    meSCV->shape_fcn(&v_shape_function_(0,0));

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_gathered_nodal_field(*coordinates_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*densityNm1_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityN_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*velocityNm1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityN_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityNp1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*Gjp_, AlgTraits::nDim_);
  dataPreReqs.add_master_element_call(SCV_VOLUME);
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
template<typename AlgTraits>
void
MomentumMassElemSuppAlg<AlgTraits>::setup()
{
  dt_ = realm_.get_time_step();
  gamma1_ = realm_.get_gamma1();
  gamma2_ = realm_.get_gamma2();
  gamma3_ = realm_.get_gamma3(); // gamma3 may be zero
}

//--------------------------------------------------------------------------
//-------- elem_execute ----------------------------------------------------
//--------------------------------------------------------------------------
template<typename AlgTraits>
void
MomentumMassElemSuppAlg<AlgTraits>::element_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity /* element */,
  ScratchViews& scratchViews
)
{
  SharedMemView<double*>& v_densityNm1 = scratchViews.get_scratch_view_1D(*densityNm1_);
  SharedMemView<double*>& v_densityN = scratchViews.get_scratch_view_1D(*densityN_);
  SharedMemView<double*>& v_densityNp1 = scratchViews.get_scratch_view_1D(*densityNp1_);
  SharedMemView<double**>& v_velocityNm1 = scratchViews.get_scratch_view_2D(*velocityNm1_);
  SharedMemView<double**>& v_velocityN = scratchViews.get_scratch_view_2D(*velocityN_);
  SharedMemView<double**>& v_velocityNp1 = scratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<double**>& v_Gpdx = scratchViews.get_scratch_view_2D(*Gjp_);

  SharedMemView<double*>& v_scv_volume = scratchViews.scv_volume;

  for (int ip=0; ip < AlgTraits::numScvIp_; ++ip) {
    const int nearestNode = ipNodeMap_[ip];

    double rhoNm1 = 0.0;
    double rhoN   = 0.0;
    double rhoNp1 = 0.0;
    for (int j=0; j < AlgTraits::nDim_; j++) {
      v_uNm1_(j) = 0.0;
      v_uN_(j) = 0.0;
      v_uNp1_(j) = 0.0;
      v_Gjp_(j) = 0.0;
    }

    for (int ic=0; ic < AlgTraits::nodesPerElement_; ++ic) {
      const double r = v_shape_function_(ip, ic);

      rhoNm1 += r * v_densityNm1(ic);
      rhoN   += r * v_densityN(ic);
      rhoNp1 += r * v_densityNp1(ic);
      for (int j=0; j < AlgTraits::nDim_; j++) {
        v_uNm1_(j) += r * v_velocityNm1(ic, j);
        v_uN_(j)   += r * v_velocityN(ic, j);
        v_uNp1_(j) += r * v_velocityNp1(ic, j);
        v_Gjp_(j)  += r * v_Gpdx(ic, j);
      }
    }

    const double scV = v_scv_volume(ip);
    const int nnNdim = nearestNode * AlgTraits::nDim_;
    // Compute RHS
    for (int j=0; j < AlgTraits::nDim_; ++j) {
      rhs[nnNdim + j] +=
        - ( gamma1_ * rhoNp1 * v_uNp1_(j) +
            gamma2_ * rhoN   * v_uN_(j) +
            gamma3_ * rhoNm1 * v_uNm1_(j)) * scV / dt_
        - v_Gjp_(j) * scV;
    }

    // Compute LHS
    const int npeNdim = AlgTraits::nodesPerElement_ * AlgTraits::nDim_;
    for (int ic=0; ic < AlgTraits::nodesPerElement_; ++ic) {
      const int icNdim = ic * AlgTraits::nDim_;
      const double r = v_shape_function_(ip, ic);
      const double lhsfac = r * gamma1_ * rhoNp1 * scV / dt_;

      for (int j=0; j<AlgTraits::nDim_; ++j) {
        const int indexNN = nnNdim + j;
        const int rowNN = indexNN * npeNdim;
        const int rNNiC_j = rowNN + icNdim + j;
        lhs[rNNiC_j] += lhsfac;
      }
    }
  }
}

INSTANTIATE_SUPPLEMENTAL_ALGORITHM(MomentumMassElemSuppAlg);

} // namespace nalu
} // namespace Sierra
