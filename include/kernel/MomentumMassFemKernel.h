/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MomentumMassFemKernel_H
#define MomentumMassFemKernel_H

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
class MomentumMassFemKernel: public Kernel
{
public:
  MomentumMassFemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    VectorFieldType*,
    ScalarFieldType*,
    ElemDataRequests&);

  virtual ~MomentumMassFemKernel();

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
  MomentumMassFemKernel() = delete;

  VectorFieldType *Gjp_{nullptr};
  VectorFieldType *coordinates_{nullptr};
  
  VectorFieldType *velocityNp1_{nullptr};
  VectorFieldType *velocityN_{nullptr};
  VectorFieldType *velocityNm1_{nullptr};
  ScalarFieldType *densityNp1_{nullptr};
  ScalarFieldType *densityN_{nullptr};
  ScalarFieldType *densityNm1_{nullptr};
 
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

#endif /* MomentumMassFemKernel_H */
