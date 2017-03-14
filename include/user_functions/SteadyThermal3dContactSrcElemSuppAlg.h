/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SteadyThermal3dContactSrcElemSuppAlg_h
#define SteadyThermal3dContactSrcElemSuppAlg_h

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
class SteadyThermal3dContactSrcElemSuppAlg : public SupplementalAlgorithm
{
public:
  SteadyThermal3dContactSrcElemSuppAlg(
    Realm &realm,
    ElemDataRequests& dataPreReqs);

  virtual ~SteadyThermal3dContactSrcElemSuppAlg() {}

  virtual void element_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    ScratchViews& scratchViews);
  
  VectorFieldType *coordinates_;

  const int *ipNodeMap_;

  const double a_;
  const double k_;
  const double pi_;

  // fixed scratch space
  Kokkos::View<double[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
  Kokkos::View<double[AlgTraits::nDim_]> v_scvCoords_{"v_scvCoords"};
};

} // namespace nalu
} // namespace Sierra

#endif
