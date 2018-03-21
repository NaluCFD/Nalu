/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SmagorinskyViscosity_h
#define SmagorinskyViscosity_h

#include <SolverAlgorithm.h>
#include <ElemDataRequests.h>
#include <FieldTypeDef.h>
#include <nalu_make_unique.h>
#include <kernel/Kernel.h>

#include <kernel/UtilityKernelManager.h>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class SmagorinskyViscosity : public NodalUtilityKernel
{
public:
  SmagorinskyViscosity(
    const stk::mesh::BulkData& bulk,
    const SolutionOptions& solnOpts,
    ElemDataRequests& requiredData,
    ElemDataRequests& outputData);

  void execute(const NodalScratchData<DoubleType>&, NodalScratchData<DoubleType>&) const;

private:
  double Cs_{0.16};
  int nDim_{3};

  stk::mesh::FieldBase* dudx_{nullptr};
  ScalarFieldType* dualVolume_{nullptr};
  ScalarFieldType* density_{nullptr};
  ScalarFieldType* turbVisc_{nullptr};

};




} // namespace nalu
} // namespace Sierra

#endif

