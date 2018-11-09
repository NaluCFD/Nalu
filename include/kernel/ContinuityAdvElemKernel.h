/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef CONTINUITYADVELEMKERNEL_H
#define CONTINUITYADVELEMKERNEL_H

#include "kernel/Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class TimeIntegrator;
class SolutionOptions;
class MasterElement;
class ElemDataRequests;

/** CMM (BDF2) for continuity equation (pressure DOF)
 */
template<typename AlgTraits>
class ContinuityAdvElemKernel: public Kernel
{
public:
  ContinuityAdvElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ElemDataRequests&);

  virtual ~ContinuityAdvElemKernel();

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
  ContinuityAdvElemKernel() = delete;

  // extract fields; nodal
  VectorFieldType *velocityRTM_{nullptr};
  VectorFieldType *Gpdx_{nullptr};
  ScalarFieldType *pressure_{nullptr};
  ScalarFieldType *densityNp1_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  double projTimeScale_{1.0};

  const bool meshMotion_;
  const bool shiftMdot_;
  const bool shiftPoisson_;
  const bool reducedSensitivities_;

  // scratch space
  AlignedViewType<DoubleType[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "view_shape_func" };

  const int* lrscv_;
};

}  // nalu
}  // sierra

#endif /* CONTINUITYADVELEMKERNEL_H */
