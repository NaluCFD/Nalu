/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ContinuityMassElemSuppAlg_h
#define ContinuityMassElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <AlgTraits.h>

// stk
#include <stk_mesh/base/Entity.hpp>

//kokkos
#include <Kokkos_Core.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;
class ElemDataRequests;
class ScratchViews;

template<typename AlgTraits>
class ContinuityMassElemSuppAlg : public SupplementalAlgorithm
{
public:

  ContinuityMassElemSuppAlg(
    Realm &realm,
    ElemDataRequests& dataPreReqs,
    const bool lumpedMass);

  virtual ~ContinuityMassElemSuppAlg() {}

  virtual void setup();

  virtual void element_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    ScratchViews& scratchViews
  );

  ScalarFieldType *densityNm1_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  VectorFieldType *coordinates_;

  double dt_;
  double gamma1_;
  double gamma2_;
  double gamma3_;
  const bool lumpedMass_;

  // master element
  const int* ipNodeMap_;

  // scratch space; geometry
  Kokkos::View<double[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_ {"view_shape_func"};
};

} // namespace nalu
} // namespace Sierra

#endif
