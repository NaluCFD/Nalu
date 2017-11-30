/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleContinuityNonConformalSolverAlgorithm_h
#define AssembleContinuityNonConformalSolverAlgorithm_h

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

class AssembleContinuityNonConformalSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleContinuityNonConformalSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *pressure,
    VectorFieldType *Gjp);
  virtual ~AssembleContinuityNonConformalSolverAlgorithm() {}

  virtual void initialize_connectivity();
  virtual void execute();

  ScalarFieldType *pressure_;
  VectorFieldType *Gjp_;
  VectorFieldType *velocity_;
  VectorFieldType *meshVelocity_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  GenericFieldType *exposedAreaVec_;
 
  const bool meshMotion_;
  
  // options that prevail over all algorithms created
  const bool useCurrentNormal_;
  const double includePstab_;
  double meshMotionFac_;

  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
