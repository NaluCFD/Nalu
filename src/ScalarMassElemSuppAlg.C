/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ScalarMassElemSuppAlg.h>
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

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// ScalarMassElemSuppAlg - CMM (BDF2/BE) for scalar equation
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
template<typename AlgTraits>
ScalarMassElemSuppAlg<AlgTraits>::ScalarMassElemSuppAlg(
  Realm &realm,
  ScalarFieldType *scalarQ,
  ElemDataRequests& dataPreReqs,
  const bool lumpedMass)
  : SupplementalAlgorithm(realm),
    scalarQNm1_(NULL),
    scalarQN_(NULL),
    scalarQNp1_(NULL),
    densityNm1_(NULL),
    densityN_(NULL),
    densityNp1_(NULL),
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
  scalarQNm1_ = realm_.number_of_states() == 2 ? &(scalarQ->field_of_state(stk::mesh::StateN)) : &(scalarQ->field_of_state(stk::mesh::StateNM1));
  scalarQN_ = &(scalarQ->field_of_state(stk::mesh::StateN));
  scalarQNp1_ = &(scalarQ->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNm1_ = realm_.number_of_states() == 2 ? &(density->field_of_state(stk::mesh::StateN)) : &(density->field_of_state(stk::mesh::StateNM1));
  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  MasterElement *meSCV = realm.get_volume_master_element(AlgTraits::topo_);

  // compute shape function
  if ( lumpedMass_ )
    meSCV->shifted_shape_fcn(&v_shape_function_(0,0));
  else
    meSCV->shape_fcn(&v_shape_function_(0,0));

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_gathered_nodal_field(*coordinates_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*scalarQNm1_, 1);
  dataPreReqs.add_gathered_nodal_field(*scalarQN_, 1);
  dataPreReqs.add_gathered_nodal_field(*scalarQNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityNm1_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityN_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_, 1);
  dataPreReqs.add_master_element_call(SCV_VOLUME);
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
template<typename AlgTraits>
void
ScalarMassElemSuppAlg<AlgTraits>::setup()
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
ScalarMassElemSuppAlg<AlgTraits>::element_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity /* element */,
  ScratchViews& scratchViews)
{
  SharedMemView<double*>& v_qNm1 = scratchViews.get_scratch_view_1D(
    *scalarQNm1_);
  SharedMemView<double*>& v_qN = scratchViews.get_scratch_view_1D(
    *scalarQN_);
  SharedMemView<double*>& v_qNp1 = scratchViews.get_scratch_view_1D(
    *scalarQNp1_);
  SharedMemView<double*>& v_rhoNm1 = scratchViews.get_scratch_view_1D(
    *densityNm1_);
  SharedMemView<double*>& v_rhoN = scratchViews.get_scratch_view_1D(
    *densityN_);
  SharedMemView<double*>& v_rhoNp1 = scratchViews.get_scratch_view_1D(
    *densityNp1_);

  SharedMemView<double*>& v_scv_volume = scratchViews.scv_volume;

  for ( int ip = 0; ip < AlgTraits::numScvIp_; ++ip ) {
      
    // nearest node to ip
    const int nearestNode = ipNodeMap_[ip];
    
    // zero out; scalar
    double qNm1Scv = 0.0;
    double qNScv = 0.0;
    double qNp1Scv = 0.0;
    double rhoNm1Scv = 0.0;
    double rhoNScv = 0.0;
    double rhoNp1Scv = 0.0;
      
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      // save off shape function
      const double r = v_shape_function_(ip,ic);

      // scalar q
      qNm1Scv += r*v_qNm1(ic);
      qNScv += r*v_qN(ic);
      qNp1Scv += r*v_qNp1(ic);

      // density
      rhoNm1Scv += r*v_rhoNm1(ic);
      rhoNScv += r*v_rhoN(ic);
      rhoNp1Scv += r*v_rhoNp1(ic);
    }

    // assemble rhs
    const double scV = v_scv_volume(ip);
    rhs[nearestNode] += 
      -(gamma1_*rhoNp1Scv*qNp1Scv + gamma2_*rhoNScv*qNScv + gamma3_*rhoNm1Scv*qNm1Scv)*scV/dt_;
    
    // manage LHS
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      // save off shape function
      const double r = v_shape_function_(ip,ic);
      const double lhsfac = r*gamma1_*rhoNp1Scv*scV/dt_;
      const int rNNiC = nearestNode*AlgTraits::nodesPerElement_+ic;
      lhs[rNNiC] += lhsfac;
    }   
  }
}
  
INSTANTIATE_SUPPLEMENTAL_ALGORITHM(ScalarMassElemSuppAlg);

} // namespace nalu
} // namespace Sierra
