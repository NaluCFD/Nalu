/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MomentumDiffFemKernel_h
#define MomentumDiffFemKernel_h

#include "kernel/Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class ElemDataRequests;
class SolutionOptions;

/** FEM diffusion kernel
 */
template<typename AlgTraits>
class MomentumDiffFemKernel: public Kernel
{
public:
  MomentumDiffFemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    VectorFieldType*,
    ScalarFieldType*,
    ElemDataRequests&);

  virtual ~MomentumDiffFemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  MomentumDiffFemKernel() = delete;

  VectorFieldType *velocityNp1_{nullptr};
  ScalarFieldType *viscosity_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  // master element
  const double includeDivU_;
  const bool shiftedGradOp_;
  
  /// Shape functions
  AlignedViewType<DoubleType[AlgTraits::numGp_]> v_ip_weight_{ "v_ip_weight" };
  AlignedViewType<DoubleType[AlgTraits::numGp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "v_shape_func" };
};

}  // nalu
}  // sierra

#endif /* MomentumDiffFemKernel_h */
