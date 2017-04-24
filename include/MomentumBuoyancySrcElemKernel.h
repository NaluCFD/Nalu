/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MOMENTUMBUOYANCYSRCELEMKERNEL_H
#define MOMENTUMBUOYANCYSRCELEMKERNEL_H

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
class ScratchViews;

/** CMM buoyancy term for momentum equation (velocity DOF)
 */
template<typename AlgTraits>
class MomentumBuoyancySrcElemKernel: public Kernel
{
public:
  MomentumBuoyancySrcElemKernel(
    const stk::mesh::BulkData&,
    SolutionOptions&,
    ElemDataRequests&);

  virtual ~MomentumBuoyancySrcElemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<double**>&,
    SharedMemView<double*>&,
    stk::mesh::Entity,
    ScratchViews&);

private:
  MomentumBuoyancySrcElemKernel() = delete;

  ScalarFieldType *densityNp1_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  double rhoRef_;
  Kokkos::View<double[AlgTraits::nDim_]> gravity_{ "view_gravity"};

  const int* ipNodeMap_;

  // scratch space
  Kokkos::View<double[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "view_shape_func" };
};

}  // nalu
}  // sierra

#endif /* MOMENTUMBUOYANCYSRCELEMKERNEL_H */
