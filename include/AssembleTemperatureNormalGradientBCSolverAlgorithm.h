/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef AssembleTemperatureNormalGradientBCSolverAlgorithm_h
#define AssembleTemperatureNormalGradientBCSolverAlgorithm_h

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

class AssembleTemperatureNormalGradientBCSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleTemperatureNormalGradientBCSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *bcScalarGrad,
    ScalarFieldType *evisc,
    ScalarFieldType *specHeat,
    bool useShifted);
  virtual ~AssembleTemperatureNormalGradientBCSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

private:

  const bool useShifted_;

  ScalarFieldType *bcScalarGrad_;
  ScalarFieldType *evisc_;
  ScalarFieldType *specHeat_;
  GenericFieldType *exposedAreaVec_;
};

}
}


#endif /* ASSEMBLESCALARELEMDIFFBCSOLVERALGORITHM_H_ */
