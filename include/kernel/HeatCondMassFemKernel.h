/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef HeatCondMassFemKernel_H
#define HeatCondMassFemKernel_H

#include "kernel/Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>
#include <vector>

namespace sierra {
namespace nalu {

class ElemDataRequests;
class SolutionOptions;
class TimeIntegrator;

/** CVFEM scalar advection/diffusion kernel
 */
template<typename AlgTraits>
class HeatCondMassFemKernel: public Kernel
{
public:
  HeatCondMassFemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ScalarFieldType*,
    ScalarFieldType*,
    ScalarFieldType*,
    ElemDataRequests&);

  virtual ~HeatCondMassFemKernel();

  /** Perform pre-timestep work for the computational kernel
   */
  virtual void setup(const TimeIntegrator&);

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  HeatCondMassFemKernel() = delete;

  const stk::mesh::BulkData* bulkData_;
  ScalarFieldType *temperatureNp1_{nullptr};
  ScalarFieldType *temperatureN_{nullptr};
  ScalarFieldType *temperatureNm1_{nullptr};
  ScalarFieldType *density_{nullptr};
  ScalarFieldType *specHeat_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  double dt_{0.0};
  double gamma1_{0.0};
  double gamma2_{0.0};
  double gamma3_{0.0};

  /// Shape functions
  AlignedViewType<DoubleType[AlgTraits::numGp_]> v_ip_weight_{ "v_ip_weight" };
  AlignedViewType<DoubleType[AlgTraits::numGp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "v_shape_func" };
};

}  // nalu
}  // sierra

#endif /* HeatCondMassFemKernel_H */
