/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumBuoyancySrcElemSuppAlg_h
#define MomentumBuoyancySrcElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <AlgTraits.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>

#include <Kokkos_Core.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;
class ElemDataRequests;
class ScratchViews;

template<typename AlgTraits>
class MomentumBuoyancySrcElemSuppAlg : public SupplementalAlgorithm
{
public:
  MomentumBuoyancySrcElemSuppAlg(
    Realm &realm,
    ElemDataRequests& dataPreReqs);

  virtual ~MomentumBuoyancySrcElemSuppAlg() {}

  virtual void element_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    ScratchViews& scratchViews);
  
  const stk::mesh::BulkData *bulkData_;

  ScalarFieldType *densityNp1_;
  VectorFieldType *coordinates_;

  double rhoRef_;
  const bool useShifted_;
  Kokkos::View<double[AlgTraits::nDim_]> gravity_{ "view_gravity"};

  const int* ipNodeMap_;

  // scratch space
  Kokkos::View<double[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "view_shape_func" };
};

} // namespace nalu
} // namespace Sierra

#endif
