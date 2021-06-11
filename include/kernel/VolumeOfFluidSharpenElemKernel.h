/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef VolumeOfFluidSharpenElemKernel_H
#define VolumeOfFluidSharpenElemKernel_H

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

template<typename AlgTraits>
class VolumeOfFluidSharpenElemKernel: public Kernel
{
public:
  VolumeOfFluidSharpenElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ScalarFieldType*,
    const double,
    ElemDataRequests&);

  virtual ~VolumeOfFluidSharpenElemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  VolumeOfFluidSharpenElemKernel() = delete;

  VectorFieldType *coordinates_{nullptr};
  VectorFieldType *velocityRTM_{nullptr};
  VectorFieldType *interfaceNormal_{nullptr};
  ScalarFieldType *vofNp1_{nullptr};

  /// Integration point to node mapping
  const double cAlpha_;
  const int* ipNodeMap_; 
};

}  // nalu
}  // sierra

#endif /* VolumeOfFluidSharpenElemKernel_H */
