/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef VOLUMETRICSRCELEMKERNEL_H
#define VOLUMETRICSRCELEMKERNEL_H

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

template<typename AlgTraits>
class VolumetricSrcElemKernel: public Kernel
{
public:
  VolumetricSrcElemKernel(
    const stk::mesh::BulkData&,
    SolutionOptions&,
    ElemDataRequests&);

  virtual ~VolumetricSrcElemKernel() {}

  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  VolumetricSrcElemKernel() = delete;

  VectorFieldType *coordinates_;

  const int *ipNodeMap_;

  const double src_;
};

}  // nalu
}  // sierra

#endif /* THERMALSRCELEMKERNEL_H */
