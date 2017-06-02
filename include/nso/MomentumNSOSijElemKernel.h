/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MOMENTUMNSOSIJELEMKERNEL_H
#define MOMENTUMNSOSIJELEMKERNEL_H

#include "Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class SolutionOptions;
class MasterElement;
class ElemDataRequests;

/** NSO for momentum equation
 *
 */
template<typename AlgTraits>
class MomentumNSOSijElemKernel : public Kernel
{
public:
  MomentumNSOSijElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    VectorFieldType*,
    ElemDataRequests&);

  virtual ~MomentumNSOSijElemKernel() {}

  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  MomentumNSOSijElemKernel() = delete;

  VectorFieldType *velocityNp1_{nullptr};
  ScalarFieldType *densityNp1_{nullptr};
  ScalarFieldType *pressure_{nullptr};
  VectorFieldType *velocityRTM_{nullptr};
  VectorFieldType *coordinates_{nullptr};
  VectorFieldType *Gjp_{nullptr};

  // master element
  const int *lrscv_;

  const double Cupw_{0.1};
  const double small_{1.0e-16};
  const double includeDivU_;
  const bool shiftedGradOp_;

  // fixed scratch space
  Kokkos::View<DoubleType[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
  Kokkos::View<DoubleType[AlgTraits::nodesPerElement_]> v_ke_{"v_ke"};
  Kokkos::View<DoubleType[AlgTraits::nDim_]> v_rhoVrtmScs_{"v_rhoVrtmScs"};
  Kokkos::View<DoubleType[AlgTraits::nDim_]> v_uNp1Scs_{"v_uNp1Scs"};
  Kokkos::View<DoubleType[AlgTraits::nDim_]> v_dpdxScs_{"v_dpdxScs"};
  Kokkos::View<DoubleType[AlgTraits::nDim_]> v_GjpScs_{"v_GjpScs"};
  Kokkos::View<DoubleType[AlgTraits::nDim_]> v_dkedxScs_{"v_dkedxScs_"};
};

}  // nalu
}  // sierra

#endif /* MOMENTUMNSOSIJELEMKERNEL_H */
