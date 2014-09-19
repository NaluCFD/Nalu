/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleRadTransWallSolverAlgorithm_h
#define AssembleRadTransWallSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;
class RadiativeTransportEquationSystem;

class AssembleRadTransWallSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleRadTransWallSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    RadiativeTransportEquationSystem *radEqSystem,
    const bool &useShifted);
  virtual ~AssembleRadTransWallSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const RadiativeTransportEquationSystem *radEqSystem_;
  const bool useShifted_;

  ScalarFieldType *intensity_;
  ScalarFieldType *bcIntensity_;
  GenericFieldType *exposedAreaVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
