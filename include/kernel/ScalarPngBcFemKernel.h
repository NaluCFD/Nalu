/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ScalarPngBcFemKernel_h
#define ScalarPngBcFemKernel_h

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

/** Wall function approach momentum equation (velocity DOF)
 */
template<typename BcAlgTraits>
class ScalarPngBcFemKernel: public Kernel
{
public:
  ScalarPngBcFemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    const std::string,
    ElemDataRequests&);

  virtual ~ScalarPngBcFemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  ScalarPngBcFemKernel() = delete;
  
  ScalarFieldType *scalarQ_{nullptr};
  
  // fixed scratch space
  AlignedViewType<DoubleType[BcAlgTraits::numFaceIp_]> v_ip_weight_{ "v_ip_weight" };
  AlignedViewType<DoubleType[BcAlgTraits::numFaceIp_][BcAlgTraits::nodesPerFace_]> vf_shape_function_{"vf_shape_function"};
};

}  // nalu
}  // sierra

#endif /* ScalarPngBcFemKernel_h */
