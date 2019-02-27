/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ScalarAdvFemKernel_h
#define ScalarAdvFemKernel_h

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
class ScalarAdvFemKernel: public Kernel
{
public:
  ScalarAdvFemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ScalarFieldType*,
    ScalarFieldType*,
    VectorFieldType*,
    ElemDataRequests&);

  virtual ~ScalarAdvFemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  ScalarAdvFemKernel() = delete;

  ScalarFieldType *scalarQ_{nullptr};
  ScalarFieldType *density_{nullptr};
  VectorFieldType *velocity_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  // master element
  const bool shiftedGradOp_;
  
  /// Shape functions
  AlignedViewType<DoubleType[AlgTraits::numGp_]> v_ip_weight_{ "v_ip_weight" };
  AlignedViewType<DoubleType[AlgTraits::numGp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "v_shape_func" };
};

}  // nalu
}  // sierra

#endif /* ScalarAdvFemKernel_h */
