/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleMeshDisplacementElemSolverAlgorithm_h
#define AssembleMeshDisplacementElemSolverAlgorithm_h

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

class AssembleMeshDisplacementElemSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleMeshDisplacementElemSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    const bool deformWrtModelCoords);
  virtual ~AssembleMeshDisplacementElemSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const bool deformWrtModelCoords_;
  VectorFieldType *meshDisplacement_;
  VectorFieldType *coordinates_;
  VectorFieldType *modelCoordinates_;
  ScalarFieldType *mu_;
  ScalarFieldType *lambda_;
};

} // namespace nalu
} // namespace Sierra

#endif
