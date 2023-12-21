/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleContinuityVofNonConformalSolverAlgorithm_h
#define AssembleContinuityVofNonConformalSolverAlgorithm_h

#include<SolutionOptions.h>
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

class AssembleContinuityVofNonConformalSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleContinuityVofNonConformalSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *pressure,
    VectorFieldType *Gjp,
    const SolutionOptions &solnOpts);
  virtual ~AssembleContinuityVofNonConformalSolverAlgorithm() {}

  virtual void initialize_connectivity();
  virtual void execute();

  ScalarFieldType *pressure_;
  VectorFieldType *Gjp_;
  VectorFieldType *velocity_;
  VectorFieldType *meshVelocity_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  ScalarFieldType *vof_;
  ScalarFieldType *interfaceCurvature_;
  ScalarFieldType *surfaceTension_;
  GenericFieldType *exposedAreaVec_;
 
  const bool meshMotion_;
  
  // options that prevail over all algorithms created
  const bool useCurrentNormal_;
  const double includePstab_;
  double meshMotionFac_;

  // local-CSF options
  const double n_;
  const double m_;
  const double c_;
  
  // added stabilization options
  std::array<double, 3> gravity_;
  double buoyancyWeight_;

  // fields to parallel communicate
  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
