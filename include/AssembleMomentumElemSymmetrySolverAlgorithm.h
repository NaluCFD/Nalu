/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleMomentumElemSymmetrySolverAlgorithm_h
#define AssembleMomentumElemSymmetrySolverAlgorithm_h

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

class AssembleMomentumElemSymmetrySolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleMomentumElemSymmetrySolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem);
  virtual ~AssembleMomentumElemSymmetrySolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const double includeDivU_;

  VectorFieldType *velocity_;
  VectorFieldType *coordinates_;
  ScalarFieldType *viscosity_;
  GenericFieldType *exposedAreaVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
