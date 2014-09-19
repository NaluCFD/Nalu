/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ASSEMBLESCALARELEMDIFFBCSOLVERALGORITHM_H_
#define ASSEMBLESCALARELEMDIFFBCSOLVERALGORITHM_H_

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

class AssembleScalarDiffBCSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleScalarDiffBCSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *bcScalarQ,
    bool use_shifted_integration);
  virtual ~AssembleScalarDiffBCSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

private:

  ScalarFieldType *bcScalarQ_;
  VectorFieldType *coordinates_;
  GenericFieldType *exposedAreaVec_;

  bool use_shifted_integration_;

};

}

}


#endif /* ASSEMBLESCALARELEMDIFFBCSOLVERALGORITHM_H_ */
