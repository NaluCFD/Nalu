/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumNSOKeElemSuppAlg_h
#define MomentumNSOKeElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <AlgTraits.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

// Kokkos
#include <Kokkos_Core.hpp>

namespace sierra{
namespace nalu{

class ElemDataRequests;
class Realm;
class ScratchViews;
class MasterElement;

template<class AlgTraits>
class MomentumNSOKeElemSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumNSOKeElemSuppAlg(
    Realm &realm,
    VectorFieldType *velocity,
    GenericFieldType *Gju,
    const double fourthFac,
    ElemDataRequests& dataPreReqs);

  virtual ~MomentumNSOKeElemSuppAlg() {}

  void element_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    ScratchViews& scratchViews);

  VectorFieldType *velocityNp1_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *pressure_;
  VectorFieldType *velocityRTM_;
  VectorFieldType *coordinates_;
  GenericFieldType *Gju_;
  VectorFieldType *Gjp_;

  // master element
  const int *lrscv_;

  const double Cupw_;
  const double small_;
  const double fourthFac_;

  // fixed scratch space
  Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
  Kokkos::View<double[AlgTraits::nodesPerElement_]> v_ke_{"v_ke"};
  Kokkos::View<double[AlgTraits::nDim_]> v_rhoVrtmScs_{"v_rhoVrtmScs"};
  Kokkos::View<double[AlgTraits::nDim_]> v_uNp1Scs_{"v_uNp1Scs"};
  Kokkos::View<double[AlgTraits::nDim_]> v_dpdxScs_{"v_dpdxScs"};
  Kokkos::View<double[AlgTraits::nDim_]> v_GjpScs_{"v_GjpScs"};
  Kokkos::View<double[AlgTraits::nDim_]> v_dkedxScs_{"v_dkedxScs_"};
};

} // namespace nalu
} // namespace Sierra

#endif
