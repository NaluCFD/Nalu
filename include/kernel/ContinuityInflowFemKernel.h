/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corp.                                           */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ContinuityInflowFemKernel_h
#define ContinuityInflowFemKernel_h

#include "FieldTypeDef.h"
#include "kernel/Kernel.h"

#include <stk_mesh/base/BulkData.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class ElemDataRequests;
class MasterElement;
class TimeIntegrator;

/** Add Int rho*uj*nj*dS
 */
template<typename BcAlgTraits>
class ContinuityInflowFemKernel: public Kernel
{
public:
  ContinuityInflowFemKernel(
    const stk::mesh::BulkData& bulkData,
    const SolutionOptions &solnOpts,
    const bool &useShifted,
    ElemDataRequests &faceDataPreReqs);

  virtual ~ContinuityInflowFemKernel();

  /** Perform pre-timestep work for the computational kernel
   */
  virtual void setup(const TimeIntegrator&);

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType **>&lhs,
    SharedMemView<DoubleType *>&rhs,
    ScratchViews<DoubleType>& scratchViews);

private:
  ContinuityInflowFemKernel() = delete;

  VectorFieldType *velocityBC_{nullptr};
  ScalarFieldType *densityBC_{nullptr};

  const bool useShifted_;
  double projTimeScale_;

  // scratch space
  AlignedViewType<DoubleType[BcAlgTraits::numFaceIp_]> v_ip_weight_{ "v_ip_weight" };
  AlignedViewType<DoubleType[BcAlgTraits::numFaceIp_][BcAlgTraits::nodesPerFace_]> v_shape_function_ { "v_shape_function" };
};

}  // nalu
}  // sierra

#endif /* ContinuityInflowFemKernel_h */
