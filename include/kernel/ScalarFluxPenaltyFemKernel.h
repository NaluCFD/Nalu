/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ScalarFluxPenaltyFemKernel_h
#define ScalarFluxPenaltyFemKernel_h

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
class ScalarFluxPenaltyFemKernel: public Kernel
{
public:
  ScalarFluxPenaltyFemKernel(
    const stk::mesh::MetaData &metaData,
    const SolutionOptions &solnOpts,
    ScalarFieldType *scalarQ,
    ScalarFieldType *bcScalarQ,
    ScalarFieldType *diffFluxCoeff,
    ElemDataRequests &faceDataPreReqs,
    ElemDataRequests &elemDataPreReqs);

  virtual ~ScalarFluxPenaltyFemKernel();

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
  ScalarFluxPenaltyFemKernel() = delete;

  ScalarFieldType *scalarQ_{nullptr};
  ScalarFieldType *bcScalarQ_{nullptr};
  ScalarFieldType *diffFluxCoeff_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  const double penaltyFac_;
  const bool shiftedGradOp_;
  MasterElement *meFEM_{nullptr};
  
  /// Shape functions
  AlignedViewType<DoubleType[BcAlgTraits::numFaceIp_]> vf_ip_weight_{ "view_face_ip_weight" };
  AlignedViewType<DoubleType[BcAlgTraits::numFaceIp_][BcAlgTraits::nodesPerFace_]> vf_shape_function_ {"view_face_shape_func"};
};

}  // nalu
}  // sierra

#endif /* ScalarFluxPenaltyFemKernel_h */
