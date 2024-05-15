/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MomentumBodyForceElemKernel_H
#define MomentumBodyForceElemKernel_H

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

template<typename AlgTraits>
class MomentumBodyForceElemKernel: public Kernel
{
public:
  MomentumBodyForceElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ElemDataRequests&);

  virtual ~MomentumBodyForceElemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  MomentumBodyForceElemKernel() = delete;
  VectorFieldType *coordinates_{nullptr};
  AlignedViewType<DoubleType[AlgTraits::nDim_]> v_body_force_{ "v_body_force" };
  const int* ipNodeMap_;
};

}  // nalu
}  // sierra

#endif /* MomentumBodyForceElemKernel_H */
