/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MomentumOpenAdvDiffElemKernel_h
#define MomentumOpenAdvDiffElemKernel_h

#include "kernel/Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class ElemDataRequests;
class MasterElement;
class PecletFunction;
class SolutionOptions;

/** Symmetry kernel for momentum equation (velocity DOF)
 */
template<typename BcAlgTraits>
class MomentumOpenAdvDiffElemKernel: public Kernel
{
public:
  MomentumOpenAdvDiffElemKernel(
    const stk::mesh::MetaData &metaData,
    const SolutionOptions &solnOpts,
    VectorFieldType *velocity,
    ScalarFieldType *viscosity,
    ElemDataRequests &faceDataPreReqs,
    ElemDataRequests &elemDataPreReqs);

  virtual ~MomentumOpenAdvDiffElemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**> &lhs,
    SharedMemView<DoubleType*> &rhs,
    ScratchViews<DoubleType> &faceScratchViews,
    ScratchViews<DoubleType> &elemScratchViews);

private:
  MomentumOpenAdvDiffElemKernel() = delete;

  GenericFieldType *Gjui{nullptr};
  ScalarFieldType *viscosity_{nullptr};
  VectorFieldType *velocityNp1_{nullptr};
  VectorFieldType *velocityRTM_{nullptr};
  VectorFieldType *coordinates_{nullptr};
  ScalarFieldType *density_{nullptr};
  GenericFieldType *exposedAreaVec_{nullptr};
  GenericFieldType *openMassFlowRate_{nullptr};
  GenericFieldType *velocityBc_{nullptr};
  
  // numerical parameters
  const double alphaUpw_;
  const double om_alphaUpw_;
  const double hoUpwind_;
  const double nfEntrain_;
  const double om_nfEntrain_;
  const double includeDivU_;
  const bool shiftedGradOp_;
  const double small_(1.0e-16);

  /// Peclet function
  PecletFunction* pecletFunction_{nullptr};

  // Integration point to node mapping
  const int* ipNodeMap_;
  const int* faceIpNodeMap_;

  /// Shape functions
  Kokkos::View<DoubleType[BcAlgTraits::numFaceIp_][BcAlgTraits::nodesPerFace_]> vf_shape_function_ {"view_face_shape_func"};
  Kokkos::View<DoubleType[BcAlgTraits::numScsIp_][BcAlgTraits::nodesPerElement_]> v_shape_function_ {"view_shape_func"};
};

}  // nalu
}  // sierra

#endif /* MomentumOpenAdvDiffElemKernel_h */
