/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ScalarNSOElemSuppAlg_h
#define ScalarNSOElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <AlgTraits.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

// Kokkos
#include <Kokkos_Core.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;
class ElemDataRequests;
class ScratchViews;

template<class AlgTraits>
class ScalarNSOElemSuppAlg : public SupplementalAlgorithm
{
public:

  ScalarNSOElemSuppAlg(
    Realm &realm,
    ScalarFieldType *scalarQ,
    VectorFieldType *Gjq,
    ScalarFieldType *diffFluxCoeff,
    const double fourthFac,
    const double altResFac,
    ElemDataRequests& dataPreReqs);

  virtual ~ScalarNSOElemSuppAlg() {}

  virtual void setup();

  virtual void element_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    ScratchViews& scratchViews);
  
  ScalarFieldType *scalarQNm1_;
  ScalarFieldType *scalarQN_;
  ScalarFieldType *scalarQNp1_;
  ScalarFieldType *densityNm1_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *velocityRTM_;
  VectorFieldType *Gjq_;
  VectorFieldType *coordinates_;

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

  // fixed space
  Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
  Kokkos::View<double[AlgTraits::nDim_]> v_dqdxScs_{"v_dukdxScs"};
  Kokkos::View<double[AlgTraits::nDim_]> v_rhoVrtmScs_{"v_rhoVrtmScs"};
};

} // namespace nalu
} // namespace Sierra

#endif
