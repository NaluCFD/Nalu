/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef VolumeOfFluidScsNoPstabUpwAdvElemKernel_H
#define VolumeOfFluidScsNoPstabUpwAdvElemKernel_H

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

/** SCS advection for VOF (upwind)
 */
template<typename AlgTraits>
class VolumeOfFluidScsNoPstabUpwAdvElemKernel: public Kernel
{
public:
  VolumeOfFluidScsNoPstabUpwAdvElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ScalarFieldType*,
    ElemDataRequests&);

  virtual ~VolumeOfFluidScsNoPstabUpwAdvElemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

  virtual DoubleType van_leer(
    const DoubleType &dqm,
    const DoubleType &dqp);

private:
  VolumeOfFluidScsNoPstabUpwAdvElemKernel() = delete;

  ScalarFieldType *vofNp1_{nullptr};
  VectorFieldType *dvofdx_{nullptr};
  VectorFieldType *velocityRTM_{nullptr};
  VectorFieldType *coordinates_{nullptr};
  
  const double hoUpwind_;
  bool useLimiter_;
  
  // Integration point to node mapping
  const int* lrscv_;
  
  const double small_{1.0e-16};

  /// Shape functions
  AlignedViewType<DoubleType[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_ {"view_shape_func"};
};

}  // nalu
}  // sierra

#endif /* VolumeOfFluidScsNoPstabUpwAdvElemKernel_H */
