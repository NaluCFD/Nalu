/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef DirichletBC_h
#define DirichletBC_h

#include <SolverAlgorithm.h>

namespace stk{
namespace mesh{
class FieldBase;
}
}

namespace sierra{
namespace nalu{

class EquationSystem;
class Realm;

class DirichletBC : public SolverAlgorithm
{
public:

  DirichletBC(
    Realm & realm,
    EquationSystem * eqSystem,
    stk::mesh::Part * part,
    stk::mesh::FieldBase * field,
    stk::mesh::FieldBase * bcValues,
    const unsigned beginPos,
    const unsigned endPos);
  
  virtual ~DirichletBC() {}

  virtual void execute();
  virtual void initialize_connectivity();

private:
  stk::mesh::FieldBase * field_;
  stk::mesh::FieldBase * bcValues_;
  const unsigned beginPos_;
  const unsigned endPos_;

private:
  // make this non-copyable
  DirichletBC(const DirichletBC & other);
  DirichletBC & operator=(const DirichletBC & other);
};

} // namespace nalu
} // namespace Sierra

#endif
