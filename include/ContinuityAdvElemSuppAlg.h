/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ContinuityAdvElemSuppAlg_h
#define ContinuityAdvElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <AlgTraits.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;
class ElemDataRequests;
class ScratchViews;

template<typename AlgTraits>
class ContinuityAdvElemSuppAlg : public SupplementalAlgorithm
{
public:

  ContinuityAdvElemSuppAlg(
    Realm &realm,
    ElemDataRequests& dataPreReqs);

  virtual ~ContinuityAdvElemSuppAlg() {}

  virtual void setup();

  virtual void element_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    ScratchViews& scratchViews
  );
  
  const stk::mesh::BulkData *bulkData_;


  // extract fields; nodal
  VectorFieldType *velocityRTM_;
  VectorFieldType *Gpdx_;
  ScalarFieldType *pressure_;
  ScalarFieldType *densityNp1_;
  VectorFieldType *coordinates_;

  double projTimeScale_;

  const bool meshMotion_;
  const bool shiftMdot_;
  const bool shiftPoisson_;
  const bool reducedSensitivities_;
  const double interpTogether_;
  const double om_interpTogether_;

  // fixed size
  Kokkos::View<double[AlgTraits::nDim_]> v_uIp_{"view_uIp"};
  Kokkos::View<double[AlgTraits::nDim_]> v_rho_uIp_{"view_rhoUIp"};
  Kokkos::View<double[AlgTraits::nDim_]> v_Gpdx_Ip_{"view_GpdxIp"};
  Kokkos::View<double[AlgTraits::nDim_]> v_dpdxIp_{"view_dpdxIp"};

  // scratch space
  Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "view_shape_func" };

  const int* lrscv_;
};

} // namespace nalu
} // namespace Sierra

#endif
