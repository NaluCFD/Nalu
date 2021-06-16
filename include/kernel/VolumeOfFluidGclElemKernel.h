/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef VolumeOfFluidGclElemKernel_H
#define VolumeOfFluidGclElemKernel_H

#include "kernel/Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class SolutionOptions;
class MasterElement;
class ElemDataRequests;

/** GCL for VOF equation (VOF DOF)
 */
template<typename AlgTraits>
class VolumeOfFluidGclElemKernel: public Kernel
{
public:
  VolumeOfFluidGclElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ScalarFieldType*,
    ElemDataRequests&,
    const bool);

  virtual ~VolumeOfFluidGclElemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  VolumeOfFluidGclElemKernel() = delete;

  ScalarFieldType *divV_{nullptr};
  ScalarFieldType *vofNp1_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  const bool lumpedMass_;

  /// Integration point to node mapping
  const int* ipNodeMap_;

  /// Shape functions
  AlignedViewType<DoubleType[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_ {"view_shape_func"};
};

}  // nalu
}  // sierra

#endif /* VolumeOfFluidGclElemKernel_H */
