/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumMassElemSuppAlg_h
#define MomentumMassElemSuppAlg_h

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
class MomentumMassElemSuppAlg : public SupplementalAlgorithm
{
public:
  MomentumMassElemSuppAlg(
    Realm &realm,
    ElemDataRequests& dataPreReqs,
    const bool lumpedMass);

  virtual ~MomentumMassElemSuppAlg() {}

  virtual void setup();

  virtual void element_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    ScratchViews& scratchViews);
  
  const stk::mesh::BulkData *bulkData_;

  VectorFieldType *velocityNm1_;
  VectorFieldType *velocityN_;
  VectorFieldType *velocityNp1_;
  ScalarFieldType *densityNm1_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  VectorFieldType *Gjp_;
  VectorFieldType *coordinates_;

  double dt_;
  double gamma1_;
  double gamma2_;
  double gamma3_;
  const bool lumpedMass_;

  // master element
  const int* ipNodeMap_;

  // scratch space
  Kokkos::View<double[AlgTraits::nDim_]> v_uNm1_ {"v_uNm1"};
  Kokkos::View<double[AlgTraits::nDim_]> v_uN_   {"v_uN"};
  Kokkos::View<double[AlgTraits::nDim_]> v_uNp1_ {"v_uNp1"};
  Kokkos::View<double[AlgTraits::nDim_]> v_Gjp_  {"v_Gjp"};

  Kokkos::View<double[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_ {"view_shape_func"};
};

} // namespace nalu
} // namespace Sierra

#endif
