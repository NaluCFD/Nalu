/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ScalarMassElemSuppAlg_h
#define ScalarMassElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <AlgTraits.h>

// stk
#include <stk_mesh/base/Entity.hpp>

// kokkos
#include <Kokkos_Core.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;
class ElemDataRequests;
class ScratchViews;

template<typename AlgTraits>
class ScalarMassElemSuppAlg : public SupplementalAlgorithm
{
public:

  ScalarMassElemSuppAlg(
    Realm &realm,
    ScalarFieldType *scalarQ,
    ElemDataRequests& dataPreReqs,
    const bool lumpedMass);

  virtual ~ScalarMassElemSuppAlg() {}

  virtual void setup();

  virtual void element_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    ScratchViews& scratchViews
  );
  
  ScalarFieldType *scalarQNm1_;
  ScalarFieldType *scalarQN_;
  ScalarFieldType *scalarQNp1_;
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
