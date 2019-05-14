/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MomentumBodyForceFemKernel_H
#define MomentumBodyForceFemKernel_H

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

/** CVFEM scalar advection/diffusion kernel
 */
template<typename AlgTraits>
class MomentumBodyForceFemKernel: public Kernel
{
public:
  MomentumBodyForceFemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ElemDataRequests&);

  virtual ~MomentumBodyForceFemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  MomentumBodyForceFemKernel() = delete;

  VectorFieldType *coordinates_{nullptr};

  /// Shape functions
  AlignedViewType<DoubleType[AlgTraits::numGp_]> v_ip_weight_{ "v_ip_weight" };
  AlignedViewType<DoubleType[AlgTraits::numGp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "v_shape_func" };
  AlignedViewType<DoubleType[AlgTraits::nDim_]> v_body_force_{ "v_body_force" };
};

}  // nalu
}  // sierra

#endif /* MomentumBodyForceFemKernel_H */
