/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ContinuityVofOpenElemKernel_h
#define ContinuityVofOpenElemKernel_h

#include "kernel/Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class SolutionOptions;
class ElemDataRequests;
class MasterElement;

/** specificed open bc (face/elem) kernel for continuity equation (pressure DOF)
 */
template<typename BcAlgTraits>
class ContinuityVofOpenElemKernel: public Kernel
{
public:
  ContinuityVofOpenElemKernel(
    const stk::mesh::MetaData &metaData,
    const SolutionOptions &solnOpts,
    ElemDataRequests &faceDataPreReqs,
    ElemDataRequests &elemDataPreReqs);

  virtual ~ContinuityVofOpenElemKernel();

  /** Perform pre-timestep work for the computational kernel
   */
  virtual void setup(const TimeIntegrator&);

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**> &lhs,
    SharedMemView<DoubleType*> &rhs,
    ScratchViews<DoubleType> &faceScratchViews,
    ScratchViews<DoubleType> &elemScratchViews,
    int elemFaceOrdinal);

private:
  ContinuityVofOpenElemKernel() = delete;

  VectorFieldType *velocityRTM_{nullptr};
  VectorFieldType *Gpdx_{nullptr};
  VectorFieldType *coordinates_{nullptr};
  ScalarFieldType *pressure_{nullptr};
  ScalarFieldType *pressureBc_{nullptr};
  ScalarFieldType *density_{nullptr};
  GenericFieldType *exposedAreaVec_{nullptr};
  GenericFieldType *dynamicPressure_{nullptr};

  double projTimeScale_{1.0};
  const double penaltyFac_{2.0};
  
  const bool shiftedGradOp_;
  const bool reducedSensitivities_;
  MasterElement *meSCS_{nullptr};

  /// Shape functions
  AlignedViewType<DoubleType[BcAlgTraits::numFaceIp_][BcAlgTraits::nodesPerFace_]> vf_shape_function_ {"view_face_shape_func"};
};

}  // nalu
}  // sierra

#endif /* ContinuityVofOpenElemKernel_h */
