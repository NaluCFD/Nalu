/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleRadTransEdgeUpwindSolverAlgorithm_h
#define AssembleRadTransEdgeUpwindSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class RadiativeTransportEquationSystem;
class Realm;

class AssembleRadTransEdgeUpwindSolverAlgorithm : public SolverAlgorithm
{
public:
  
  AssembleRadTransEdgeUpwindSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    RadiativeTransportEquationSystem *radEqSystem);
  virtual ~AssembleRadTransEdgeUpwindSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const RadiativeTransportEquationSystem *radEqSystem_;

  ScalarFieldType *intensity_;
  VectorFieldType *edgeAreaVec_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
