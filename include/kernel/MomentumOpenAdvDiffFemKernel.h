/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MomentumOpenAdvDiffFemKernel_h
#define MomentumOpenAdvDiffFemKernel_h

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
class MomentumOpenAdvDiffFemKernel: public Kernel
{
public:
  MomentumOpenAdvDiffFemKernel(
    const stk::mesh::MetaData &metaData,
    const SolutionOptions &solnOpts,
    VectorFieldType *,
    ScalarFieldType *,
    ElemDataRequests &faceDataPreReqs,
    ElemDataRequests &elemDataPreReqs);

  virtual ~MomentumOpenAdvDiffFemKernel();

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
  MomentumOpenAdvDiffFemKernel() = delete;

  VectorFieldType *velocityNp1_{nullptr};
  VectorFieldType *vrtmL_{nullptr};
  VectorFieldType *velocityBc_{nullptr};
  VectorFieldType *meshVelocity_{nullptr};
  VectorFieldType *GjpL_{nullptr};
  ScalarFieldType *pressure_{nullptr};
  ScalarFieldType *pressureBc_{nullptr};
  ScalarFieldType *density_{nullptr};
  ScalarFieldType *viscosity_{nullptr};
  GenericFieldType *dynamicPressure_{nullptr};

  const double includeDivU_;
  const double meshVelocityCorrection_;
  const bool shiftedGradOp_;
  double projTimeScale_;
  const double penaltyFac_;

  MasterElement *meFEM_{nullptr};

  // fixed scratch space
  AlignedViewType<DoubleType[BcAlgTraits::numFaceIp_]> vf_ip_weight_{ "view_face_ip_weight" };
  AlignedViewType<DoubleType[BcAlgTraits::numFaceIp_][BcAlgTraits::nodesPerFace_]> vf_shape_function_ {"view_face_shape_func"};
};

}  // nalu
}  // sierra

#endif /* MomentumOpenAdvDiffFemKernel_h */
