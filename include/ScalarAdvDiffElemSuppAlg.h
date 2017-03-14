/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ScalarAdvDiffElemSuppAlg_h
#define ScalarAdvDiffElemSuppAlg_h

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
class ScalarAdvDiffElemSuppAlg : public SupplementalAlgorithm
{
public:
  ScalarAdvDiffElemSuppAlg(
    Realm &realm,
    ScalarFieldType *scalarQ,
    ScalarFieldType *diffFluxCoeff,
    ElemDataRequests& dataPreReqs);

  virtual ~ScalarAdvDiffElemSuppAlg() {}

  void element_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    ScratchViews& scratchViews);
  
  ScalarFieldType *scalarQ_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *coordinates_;
  GenericFieldType *massFlowRate_;

  // master element
  const int *lrscv_;

  // fixed scratch space
  Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
};

} // namespace nalu
} // namespace Sierra

#endif
