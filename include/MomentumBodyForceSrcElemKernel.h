/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MOMENTUMBODYFORCESRCELEMKERNEL_H
#define MOMENTUMBODYFORCESRCELEMKERNEL_H

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

/** CMM BodyForce term for momentum equation (velocity DOF)
 */
template<typename AlgTraits>
class MomentumBodyForceSrcElemKernel: public Kernel
{
public:
  MomentumBodyForceSrcElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ElemDataRequests&);

  virtual ~MomentumBodyForceSrcElemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<double**>&,
    SharedMemView<double*>&,
    stk::mesh::Entity,
    ScratchViews&);

private:
  MomentumBodyForceSrcElemKernel() = delete;

  VectorFieldType *coordinates_{nullptr};

  Kokkos::View<double[AlgTraits::nDim_]> v_body_force_{ "v_body_force"};

  const int* ipNodeMap_;
};

}  // nalu
}  // sierra

#endif /* MOMENTUMBODYFORCESRCELEMKERNEL_H */
