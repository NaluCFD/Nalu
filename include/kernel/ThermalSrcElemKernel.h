/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef THERMALSRCELEMKERNEL_H
#define THERMALSRCELEMKERNEL_H

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
class ThermalSrcElemKernel: public Kernel
{
public:
  ThermalSrcElemKernel(
    const stk::mesh::BulkData&,
    SolutionOptions&,
    ElemDataRequests&);

  virtual ~ThermalSrcElemKernel() {}

  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  ThermalSrcElemKernel() = delete;

  VectorFieldType *coordinates_;

  const int *ipNodeMap_;

  const double src_;
};

}  // nalu
}  // sierra

#endif /* THERMALSRCELEMKERNEL_H */
