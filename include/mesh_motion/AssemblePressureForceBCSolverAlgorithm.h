/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef AssemblePressureForceBCSolverAlgorithm_h
#define AssemblePressureForceBCSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class LinearSystem;
class Realm;

class AssemblePressureForceBCSolverAlgorithm : public SolverAlgorithm
{
public:

  AssemblePressureForceBCSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *bcScalarQ,
    bool use_shifted_integration);
  virtual ~AssemblePressureForceBCSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

private:

  ScalarFieldType *bcScalarQ_;
  VectorFieldType *coordinates_;
  GenericFieldType *exposedAreaVec_;

  bool use_shifted_integration_;

};

} // namespace nalu
} // namespace Sierra

#endif
