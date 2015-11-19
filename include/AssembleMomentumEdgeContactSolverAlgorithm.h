/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleMomentumEdgeContactSolverAlgorithm_h
#define AssembleMomentumEdgeContactSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
class FieldBase;
}
}

namespace sierra{
namespace nalu{

class Realm;
class PecletFunction;

class AssembleMomentumEdgeContactSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleMomentumEdgeContactSolverAlgorithm(
      Realm &realm,
      stk::mesh::Part *part,
      EquationSystem *eqSystem);
  virtual ~AssembleMomentumEdgeContactSolverAlgorithm();
  virtual void initialize_connectivity();
  virtual void execute();

  double van_leer(
      const double &dqm,
      const double &dqp,
      const double &small);

  const bool meshMotion_;
  const double includeDivU_;
  
  VectorFieldType *meshVelocity_;
  VectorFieldType *velocity_;
  VectorFieldType *coordinates_;
  GenericFieldType *dudx_;
  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;
  ScalarFieldType *haloMdot_;

  // peclet function specifics
  PecletFunction * pecletFunction_;

  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
