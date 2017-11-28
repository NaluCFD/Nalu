/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef AssembleScalarFluxBCSolverAlgorithm_h
#define AssembleScalarFluxBCSolverAlgorithm_h

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

class AssembleScalarFluxBCSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleScalarFluxBCSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *bcScalarQ,
    bool useShifted);
  virtual ~AssembleScalarFluxBCSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

private:

  const bool useShifted_;

  ScalarFieldType *bcScalarQ_;
  GenericFieldType *exposedAreaVec_;
};

}
}


#endif /* ASSEMBLESCALARELEMDIFFBCSOLVERALGORITHM_H_ */
