/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumAdvDiffElemSuppAlg_h
#define MomentumAdvDiffElemSuppAlg_h

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
class MomentumAdvDiffElemSuppAlg : public SupplementalAlgorithm
{
public:
 
  MomentumAdvDiffElemSuppAlg(
    Realm &realm,
    VectorFieldType *velocity,
    ScalarFieldType *viscosity,
    ElemDataRequests& dataPreReqs);

  virtual ~MomentumAdvDiffElemSuppAlg() {}

  void element_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    ScratchViews& scratchViews);

  VectorFieldType *velocityNp1_;
  VectorFieldType *coordinates_;
  ScalarFieldType *viscosity_;
  GenericFieldType *massFlowRate_;

  // master element
  const int *lrscv_;

  const double includeDivU_;

  // fixed scratch space
  Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
  Kokkos::View<double[AlgTraits::nDim_]> v_uIp_{"v_uIp"};
};

} // namespace nalu
} // namespace Sierra

#endif
