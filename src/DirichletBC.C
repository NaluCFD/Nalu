/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <DirichletBC.h>
#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <Realm.h>

namespace sierra{
namespace nalu{

DirichletBC::DirichletBC(
  Realm & realm,
  EquationSystem * eqSystem,
  stk::mesh::Part * part,
  stk::mesh::FieldBase * field,
  stk::mesh::FieldBase * bcValues,
  const unsigned beginPos,
  const unsigned endPos)
  : SolverAlgorithm(realm, part, eqSystem),
    field_(field),
    bcValues_(bcValues),
    beginPos_(beginPos),
    endPos_(endPos)
{}

void
DirichletBC::initialize_connectivity()
{
  eqSystem_->linsys_->buildDirichletNodeGraph(partVec_);
}

void
DirichletBC::execute()
{

  eqSystem_->linsys_->applyDirichletBCs(
    field_,
    bcValues_,
    partVec_,
    beginPos_,
    endPos_);
}

} // namespace nalu
} // namespace Sierra
