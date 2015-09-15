/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssemblePNGElemSolverAlgorithm_h
#define AssemblePNGElemSolverAlgorithm_h

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

class AssemblePNGElemSolverAlgorithm : public SolverAlgorithm
{
public:

  AssemblePNGElemSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    std::string independentDofName,
    std::string dofName);
  virtual ~AssemblePNGElemSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  ScalarFieldType *scalarQ_;
  VectorFieldType *dqdx_;
  VectorFieldType *coordinates_;
};

} // namespace nalu
} // namespace Sierra

#endif
