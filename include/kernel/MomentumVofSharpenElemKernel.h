/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MomentumVofSharpenElemKernel_H
#define MomentumVofSharpenElemKernel_H

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

/** Tanh sharpening kernel
 */
template<typename AlgTraits>
class MomentumVofSharpenElemKernel: public Kernel
{
public:
  MomentumVofSharpenElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    VectorFieldType*,
    ElemDataRequests&);

  virtual ~MomentumVofSharpenElemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  MomentumVofSharpenElemKernel() = delete;

  VectorFieldType *coordinates_{nullptr};
  VectorFieldType *velocityNp1_{nullptr};
  VectorFieldType *velocityRTM_{nullptr};
  ScalarFieldType *vof_{nullptr};
  VectorFieldType *interfaceNormal_{nullptr};

  // integration point to node mapping
  const int* lrscv_;

  // constants
  const double cAlpha_;
  const double densDiff_;

  // fixed scratch space
  AlignedViewType<DoubleType[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
};

}  // nalu
}  // sierra

#endif /* MomentumVofSharpenElemKernel_H */
