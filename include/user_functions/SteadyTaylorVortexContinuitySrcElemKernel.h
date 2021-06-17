/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SteadyTaylorVortexContinuitySrcElemKernel_H
#define SteadyTaylorVortexContinuitySrcElemKernel_H

#include "kernel/Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class TimeIntegrator;
class SolutionOptions;
class MasterElement;
class ElemDataRequests;

template<typename AlgTraits>
class SteadyTaylorVortexContinuitySrcElemKernel : public Kernel
{
public:
  SteadyTaylorVortexContinuitySrcElemKernel(
    const stk::mesh::BulkData&,
    SolutionOptions&,
    ElemDataRequests&);

  virtual ~SteadyTaylorVortexContinuitySrcElemKernel() {}

  /** Perform pre-timestep work for the computational kernel
   */
  virtual void setup(const TimeIntegrator&);

  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  SteadyTaylorVortexContinuitySrcElemKernel() = delete;

  VectorFieldType *coordinates_;

  const int *ipNodeMap_;
  const double rhoP_;
  const double rhoS_;
  const double unot_;
  const double vnot_;
  const double znot_;
  const double pnot_;
  const double a_;
  const double amf_;
  const double Sc_;  
  const double pi_;
  double projTimeScale_;

  // fixed scratch space
  AlignedViewType<DoubleType[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
};

}  // nalu
}  // sierra

#endif /* SteadyTaylorVortexContinuitySrcElemKernel_H */
