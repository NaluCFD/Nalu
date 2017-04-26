/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SCALARNSOELEMKERNEL_H
#define SCALARNSOELEMKERNEL_H

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
class ScratchViews;

/** NSO for momentum equation
 *
 */
template<typename AlgTraits>
class ScalarNSOElemKernel : public Kernel
{
public:
  ScalarNSOElemKernel(
    const stk::mesh::BulkData&,
    SolutionOptions&,
    ScalarFieldType*,
    VectorFieldType*,
    ScalarFieldType*,
    const double,
    const double,
    ElemDataRequests&);

  virtual ~ScalarNSOElemKernel() {}

  virtual void setup(const TimeIntegrator&);

  virtual void execute(
    SharedMemView<double**>&,
    SharedMemView<double*>&,
    stk::mesh::Entity,
    ScratchViews&);

private:
  ScalarNSOElemKernel() = delete;

  ScalarFieldType *scalarQNm1_{nullptr};
  ScalarFieldType *scalarQN_{nullptr};
  ScalarFieldType *scalarQNp1_{nullptr};
  ScalarFieldType *densityNm1_{nullptr};
  ScalarFieldType *densityN_{nullptr};
  ScalarFieldType *densityNp1_{nullptr};
  ScalarFieldType *diffFluxCoeff_{nullptr};
  VectorFieldType *velocityRTM_{nullptr};
  VectorFieldType *Gjq_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  const int *lrscv_;

  double dt_{0.0};
  double gamma1_{0.0};
  double gamma2_{0.0};
  double gamma3_{0.0};
  const double Cupw_{0.1};
  const double fourthFac_;
  const double altResFac_;
  const double om_altResFac_;
  const double nonConservedForm_{0.0};
  const double small_{1.0e-16};

  // fixed scratch space
  Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
  Kokkos::View<double[AlgTraits::nDim_]> v_dqdxScs_{"v_dqdxScs"};
  Kokkos::View<double[AlgTraits::nDim_]> v_rhoVrtmScs_{"v_rhoVrtmScs"};
};

}  // nalu
}  // sierra

#endif /* SCALARNSOELEMKERNEL_H */
