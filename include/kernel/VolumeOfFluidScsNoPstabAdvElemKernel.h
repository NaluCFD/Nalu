/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef VolumeOfFluidScsNoPstabAdvElemKernel_H
#define VolumeOfFluidScsNoPstabAdvElemKernel_H

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

/** SCS advection for VOF
 */
template<typename AlgTraits>
class VolumeOfFluidScsNoPstabAdvElemKernel: public Kernel
{
public:
  VolumeOfFluidScsNoPstabAdvElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ScalarFieldType*,
    ElemDataRequests&);

  virtual ~VolumeOfFluidScsNoPstabAdvElemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  VolumeOfFluidScsNoPstabAdvElemKernel() = delete;

  ScalarFieldType *vofNp1_{nullptr};
  VectorFieldType *coordinates_{nullptr};
  GenericFieldType *velocityRTM_{nullptr};

  // Integration point to node mapping
  const int* lrscv_;

  /// Shape functions
  AlignedViewType<DoubleType[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_ {"view_shape_func"};
};

}  // nalu
}  // sierra

#endif /* VolumeOfFluidScsNoPstabAdvElemKernel_H */
