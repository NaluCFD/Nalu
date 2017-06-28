/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef STEADYTHERMAL3DCONTACTSRCELEMKERNEL_H
#define STEADYTHERMAL3DCONTACTSRCELEMKERNEL_H

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

template<typename AlgTraits>
class SteadyThermal3dContactSrcElemKernel: public Kernel
{
public:
  SteadyThermal3dContactSrcElemKernel(
    const stk::mesh::BulkData&,
    SolutionOptions&,
    ElemDataRequests&);

  virtual ~SteadyThermal3dContactSrcElemKernel() {}

  virtual void execute(
    SharedMemView<double**>&,
    SharedMemView<double*>&,
    ScratchViews&);

private:
  SteadyThermal3dContactSrcElemKernel() = delete;

  VectorFieldType *coordinates_;

  const int *ipNodeMap_;

  const double a_;
  const double k_;
  const double pi_;

  // fixed scratch space
  Kokkos::View<double[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
};

}  // nalu
}  // sierra

#endif /* STEADYTHERMAL3DCONTACTSRCELEMKERNEL_H */
