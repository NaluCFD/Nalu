/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef VolumeOfFluidOpenAdvElemKernel_h
#define VolumeOfFluidOpenAdvElemKernel_h

#include "master_element/MasterElement.h"

// scratch space
#include "ScratchViews.h"

#include "kernel/Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class ElemDataRequests;
class EquationSystem;
class MasterElement;
class SolutionOptions;

/** open advection for vof
 */
template<typename BcAlgTraits>
class VolumeOfFluidOpenAdvElemKernel: public Kernel
{
public:
  VolumeOfFluidOpenAdvElemKernel(
    const stk::mesh::MetaData &metaData,
    const SolutionOptions &solnOpts,
    EquationSystem* eqSystem,  
    ScalarFieldType *scalarQ,
    ScalarFieldType *bcScalarQ,
    ElemDataRequests &faceDataPreReqs,
    ElemDataRequests &elemDataPreReqs);

  virtual ~VolumeOfFluidOpenAdvElemKernel();

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
  VolumeOfFluidOpenAdvElemKernel() = delete;

  ScalarFieldType *scalarQ_{nullptr};
  ScalarFieldType *bcScalarQ_{nullptr};
  VectorFieldType *coordinates_{nullptr};
  GenericFieldType *openVolumeFlowRate_{nullptr};
  
  // Integration point to node mapping
  MasterElement *meSCS_{nullptr};

  /// Shape functions
  AlignedViewType<DoubleType[BcAlgTraits::numFaceIp_][BcAlgTraits::nodesPerFace_]> vf_adv_shape_function_ {"vf_adv_shape_function"};
};

}  // nalu
}  // sierra

#endif /* VolumeOfFluidOpenAdvElemKernel_h */
