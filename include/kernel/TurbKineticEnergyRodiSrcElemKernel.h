/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef TurbKineticEnergyRodiSrcElemKernel_H
#define TurbKineticEnergyRodiSrcElemKernel_H

#include "kernel/Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class SolutionOptions;
class MasterElement;
class ElemDataRequests;

/** Add Rodi source term for kernel-based algorithm approach
 */
template<typename AlgTraits>
class TurbKineticEnergyRodiSrcElemKernel: public Kernel
{
public:
  TurbKineticEnergyRodiSrcElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ElemDataRequests&);

  virtual ~TurbKineticEnergyRodiSrcElemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  TurbKineticEnergyRodiSrcElemKernel() = delete;

  VectorFieldType *coordinates_{nullptr};
  VectorFieldType *dhdx_{nullptr};
  ScalarFieldType *specificHeat_{nullptr};
  ScalarFieldType *tvisc_{nullptr};

  const double beta_;
  const double turbPr_;
  
  // Integration point to node mapping
  const int* ipNodeMap_;

  AlignedViewType<DoubleType[AlgTraits::nDim_]> gravity_{ "v_gravity"};
};

}  // nalu
}  // sierra

#endif /* TurbKineticEnergyRodiSrcElemKernel_H */
