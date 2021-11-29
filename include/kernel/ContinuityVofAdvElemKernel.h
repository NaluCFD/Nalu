/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ContinuityVofAdvElemKernel_H
#define ContinuityVofAdvElemKernel_H

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

/** CMM (BDF2) for continuity equation (pressure DOF)
 */
template<typename AlgTraits>
class ContinuityVofAdvElemKernel: public Kernel
{
public:
  ContinuityVofAdvElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ElemDataRequests&);

  virtual ~ContinuityVofAdvElemKernel();

  /** Perform pre-timestep work for the computational kernel
   */
  virtual void setup(const TimeIntegrator&);

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  ContinuityVofAdvElemKernel() = delete;

  // extract fields; nodal
  VectorFieldType *velocityRTM_{nullptr};
  VectorFieldType *Gpdx_{nullptr};
  ScalarFieldType *pressure_{nullptr};
  ScalarFieldType *densityNp1_{nullptr};
  ScalarFieldType *interfaceCurvature_{nullptr};
  ScalarFieldType *surfaceTension_{nullptr};
  ScalarFieldType *vof_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  double projTimeScale_{1.0};

  const bool meshMotion_;
  const bool shiftMdot_;
  const bool shiftPoisson_;
  const bool reducedSensitivities_;

  AlignedViewType<DoubleType[AlgTraits::nDim_]> gravity_{ "v_gravity"};

  // scratch space
  AlignedViewType<DoubleType[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "view_shape_func" };

  const int* lrscv_;

  DoubleType buoyancyWeight_{0.0};

};

}  // nalu
}  // sierra

#endif /* ContinuityVofAdvElemKernel_H */
