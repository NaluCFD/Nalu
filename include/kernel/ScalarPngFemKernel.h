/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ScalarPngFemKernel_H
#define ScalarPngFemKernel_H

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

/** FEM PNG kernel
 */
template<typename AlgTraits>
class ScalarPngFemKernel: public Kernel
{
public:
  ScalarPngFemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    std::string,
    std::string,
    ElemDataRequests&);

  virtual ~ScalarPngFemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  ScalarPngFemKernel() = delete;

  ScalarFieldType *scalarQ_{nullptr};
  VectorFieldType *Gjq_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  /// Shape functions
  AlignedViewType<DoubleType[AlgTraits::numGp_]> v_ip_weight_{ "v_ip_weight" };
  AlignedViewType<DoubleType[AlgTraits::numGp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "v_shape_func" };
};

}  // nalu
}  // sierra

#endif /* ScalarPngFemKernel_H */
