/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef CONTINUITYADVFEMKERNEL_H
#define CONTINUITYADVFEMKERNEL_H

#include "Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class ElemDataRequests;
class SolutionOptions;
class TimeIntegrator;

/** CVFEM scalar advection/diffusion kernel
 */
template<typename AlgTraits>
class ContinuityAdvFemKernel: public Kernel
{
public:
  ContinuityAdvFemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ElemDataRequests&);

  virtual ~ContinuityAdvFemKernel();

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
  ContinuityAdvFemKernel() = delete;

  VectorFieldType *velocityRTM_{nullptr};
  VectorFieldType *Gjp_{nullptr};
  ScalarFieldType *pressure_{nullptr};
  ScalarFieldType *densityNp1_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  // master element
  const bool meshMotion_;
  const bool shiftMdot_;
  const bool shiftedGradOp_;
  const bool reducedSensitivities_;
  double projTimeScale_;
  
  /// Shape functions
  AlignedViewType<DoubleType[AlgTraits::numGp_]> v_ip_weight_{ "v_ip_weight" };
  AlignedViewType<DoubleType[AlgTraits::numGp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "v_shape_func" };
};

}  // nalu
}  // sierra

#endif /* CONTINUITYADVFEMKERNEL_H */
