/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SCALARDIFFELEMKERNEL_H
#define SCALARDIFFELEMKERNEL_H

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

/** CVFEM scalar advection/diffusion kernel
 */
template<typename AlgTraits>
class ScalarDiffElemKernel: public Kernel
{
public:
  ScalarDiffElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ScalarFieldType*,
    ScalarFieldType*,
    ElemDataRequests&);

  virtual ~ScalarDiffElemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<double**>&,
    SharedMemView<double*>&,
    stk::mesh::Entity,
    ScratchViews&);

private:
  ScalarDiffElemKernel() = delete;

  ScalarFieldType *scalarQ_{nullptr};
  ScalarFieldType *diffFluxCoeff_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  /// Left right node indicators
  const int* lrscv_;

  /// Shape functions
  Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "v_shape_func" };
};

}  // nalu
}  // sierra

#endif /* SCALARDIFFELEMKERNEL_H */
