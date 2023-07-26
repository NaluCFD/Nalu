/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MomentumVofCapillaryElemKernel_H
#define MomentumVofCapillaryElemKernel_H

#include "kernel/Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class SolutionOptions;
class TimeIntegrator;
class MasterElement;
class ElemDataRequests;

/** capillary stabilization kernel for the momentum equation (velocity DOF) 
 */
template<typename AlgTraits>
class MomentumVofCapillaryElemKernel: public Kernel
{
public:
  MomentumVofCapillaryElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    VectorFieldType*,
    ElemDataRequests&);

  virtual ~MomentumVofCapillaryElemKernel();

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
  MomentumVofCapillaryElemKernel() = delete;

  VectorFieldType *velocityNp1_{nullptr};
  VectorFieldType *coordinates_{nullptr};
  VectorFieldType *interfaceNormal_{nullptr};
  ScalarFieldType *vof_{nullptr};
  ScalarFieldType *surfaceTension_{nullptr};

  double dt_{0.0};

  const int* lrscv_;

  const bool shiftedGradOp_;

  // fixed scratch space
  AlignedViewType<DoubleType[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
};

}  // nalu
}  // sierra

#endif /* MOMENTUMADVDIFFELEMKERNEL_H */
