/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssemblePNGBoundarySolverAlgorithm_h
#define AssemblePNGBoundarySolverAlgorithm_h

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

class AssemblePNGBoundarySolverAlgorithm : public SolverAlgorithm
{
public:

  AssemblePNGBoundarySolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    std::string independentDofName);
  virtual ~AssemblePNGBoundarySolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  ScalarFieldType *scalarQ_;
  GenericFieldType *exposedAreaVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
