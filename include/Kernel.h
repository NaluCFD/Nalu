/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef KERNEL_H
#define KERNEL_H

#include "KokkosInterface.h"
#include "SimdInterface.h"
#include "ScratchViews.h"

#include <stk_mesh/base/Entity.hpp>

namespace sierra {
namespace nalu {

class TimeIntegrator;
class SolutionOptions;

/** Base class for computational kernels in Nalu
 *
 * A kernel represents an atomic unit of computation applied over a given set of
 * nodes, elements, etc. using STK and Kokkos looping constructs.
 */
class Kernel
{
public:
  Kernel() = default;

  virtual ~Kernel() {}

  /** Perform pre-timestep work for the computational kernel
   */
  virtual void setup(const TimeIntegrator&) {}

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&)
  {}
};

}  // nalu
}  // sierra

#endif /* KERNEL_H */
