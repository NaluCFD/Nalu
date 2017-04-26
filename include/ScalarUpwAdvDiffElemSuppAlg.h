/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ScalarUpwAdvDiffElemSuppAlg_h
#define ScalarUpwAdvDiffElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <AlgTraits.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

// Kokkos
#include <Kokkos_Core.hpp>

namespace sierra{
namespace nalu{

class EquationSystem;
class ElemDataRequests;
class Realm;
class PecletFunction;
class ScratchViews;
class MasterElement;

template<class AlgTraits>
class ScalarUpwAdvDiffElemSuppAlg : public SupplementalAlgorithm
{
public:
  ScalarUpwAdvDiffElemSuppAlg(
    Realm &realm,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    VectorFieldType *Gjq,
    ScalarFieldType *diffFluxCoeff,
    ElemDataRequests& dataPreReqs);

  virtual ~ScalarUpwAdvDiffElemSuppAlg();

  virtual void element_execute(
    SharedMemView<double **>& lhs,
    SharedMemView<double *>& rhs,
    stk::mesh::Entity element,
    ScratchViews& scratchViews);

  virtual void setup();

  virtual double van_leer(
    const double &dqm,
    const double &dqp);

  ScalarFieldType *scalarQ_;
  VectorFieldType *Gjq_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *velocityRTM_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  GenericFieldType *massFlowRate_;

  // master element
  const int *lrscv_;

  // upwind options
  const std::string dofName_;
  double alpha_;
  double alphaUpw_;
  double hoUpwind_;
  bool useLimiter_;
  double om_alpha_;
  double om_alphaUpw_;
  const double small_;

  // peclet function specifics
  PecletFunction * pecletFunction_;

  // fixed scratch space
  Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
  Kokkos::View<double[AlgTraits::nDim_]> v_coordIp_{"v_coordIp"};
};

} // namespace nalu
} // namespace Sierra

#endif
