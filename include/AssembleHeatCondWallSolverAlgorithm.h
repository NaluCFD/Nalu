/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleHeatCondWallSolverAlgorithm_h
#define AssembleHeatCondWallSolverAlgorithm_h

#include <SolverAlgorithm.h>
#include <FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleHeatCondWallSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleHeatCondWallSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *referenceTemp,
    ScalarFieldType *couplingParameter,
    ScalarFieldType *normalHeatFlux,
    bool useShifted = false);
  virtual ~AssembleHeatCondWallSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const bool useShifted_;

  GenericFieldType *exposedAreaVec_;
  ScalarFieldType *referenceTemp_;
  ScalarFieldType *couplingParameter_;
  ScalarFieldType *normalHeatFlux_;
  ScalarFieldType *temperature_;
};

} // namespace nalu
} // namespace Sierra

#endif
