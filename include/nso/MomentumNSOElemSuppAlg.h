/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumNSOElemSuppAlg_h
#define MomentumNSOElemSuppAlg_h

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
class MomentumNSOElemSuppAlg : public SupplementalAlgorithm
{
public:
 
  MomentumNSOElemSuppAlg(
    Realm &realm,
    VectorFieldType *velocity,
    GenericFieldType *Gju,
    ScalarFieldType *viscosity,
    const double fourthFac,
    const double altResFac,
    ElemDataRequests& dataPreReqs);

  virtual ~MomentumNSOElemSuppAlg() {}

  void setup();

  void element_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    ScratchViews& scratchViews);

  VectorFieldType *velocityNm1_;
  VectorFieldType *velocityN_;
  VectorFieldType *velocityNp1_;
  ScalarFieldType *densityNm1_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *pressure_;
  VectorFieldType *velocityRTM_;
  VectorFieldType *coordinates_;
  ScalarFieldType *viscosity_;
  GenericFieldType *Gju_;

  // master element
  const int *lrscv_;

  double dt_;
  double gamma1_;
  double gamma2_;
  double gamma3_;
  const double Cupw_;
  const double small_;
  const double fourthFac_;
  const double altResFac_;
  const double om_altResFac_;
  const double nonConservedForm_;
  const double includeDivU_;

  // fixed scratch space
  Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
  Kokkos::View<double[AlgTraits::nDim_]> v_dukdxScs_{"v_dukdxScs"};
  Kokkos::View<double[AlgTraits::nDim_]> v_rhoVrtmScs_{"v_rhoVrtmScs"};
  Kokkos::View<double[AlgTraits::nDim_]> v_dpdxScs_{"v_dpdxScs"};
  Kokkos::View<double[AlgTraits::nDim_][AlgTraits::nDim_]> v_kd_{"v_kd"};
};

} // namespace nalu
} // namespace Sierra

#endif
