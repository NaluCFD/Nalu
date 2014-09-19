/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleHeatCondIrradWallSolverAlgorithm_h
#define AssembleHeatCondIrradWallSolverAlgorithm_h

#include <SolverAlgorithm.h>
#include <FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleHeatCondIrradWallSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleHeatCondIrradWallSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    bool useShifted = false);
  virtual ~AssembleHeatCondIrradWallSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const bool useShifted_;

  GenericFieldType *exposedAreaVec_;
  ScalarFieldType *temperature_;
  ScalarFieldType *irradiation_;
  ScalarFieldType *emissivity_;
};

} // namespace nalu
} // namespace Sierra

#endif
