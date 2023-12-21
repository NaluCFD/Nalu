/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleVofNonConformalSolverAlgorithm_h
#define AssembleVofNonConformalSolverAlgorithm_h

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

class AssembleVofNonConformalSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleVofNonConformalSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *vof);
  virtual ~AssembleVofNonConformalSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  ScalarFieldType *vof_;
  VectorFieldType *coordinates_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *ncVolumeFlowRate_;

  // options that prevail over all algorithms created
  const double eta_;
  const bool useCurrentNormal_;

  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
