/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef TurbKineticEnergyKsgsBuoyancySrcElemKernel_H
#define TurbKineticEnergyKsgsBuoyancySrcElemKernel_H

#include "Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class SolutionOptions;
class MasterElement;
class ElemDataRequests;

/** Add Ksgs source term for kernel-based algorithm approach
 */
template<typename AlgTraits>
class TurbKineticEnergyKsgsBuoyancySrcElemKernel: public Kernel
{
public:
  TurbKineticEnergyKsgsBuoyancySrcElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ElemDataRequests&);

  virtual ~TurbKineticEnergyKsgsBuoyancySrcElemKernel();
  
  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  TurbKineticEnergyKsgsBuoyancySrcElemKernel() = delete;
  
  VectorFieldType *coordinates_{nullptr};
  ScalarFieldType *tkeNp1_{nullptr};
  ScalarFieldType *densityNp1_{nullptr};
  ScalarFieldType *dualNodalVolume_{nullptr};
  
  const double CbTwo_;

  // Integration point to node mapping
  const int* ipNodeMap_;
  
  // drho/dx X gravity
  virtual DoubleType cross_product_magnitude();

  // fixed scratch space; non-view
  double gravity_[3];

  // view
  AlignedViewType<DoubleType[AlgTraits::nDim_]> w_drhodx_{"w_drhodx"};
  AlignedViewType<DoubleType[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
};
 
}  // nalu
}  // sierra

#endif /* TurbKineticEnergyKsgsBuoyancySrcElemKernel_H */
