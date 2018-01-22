/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleElemSolverAlgorithm_h
#define AssembleElemSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<ElemDataRequests.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
class Topology;
}
}

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class AssembleElemSolverAlgorithm : public SolverAlgorithm
{
public:
  AssembleElemSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    stk::mesh::EntityRank entityRank,
    unsigned nodesPerEntity,
    bool interleaveMeViews = true);
  virtual ~AssembleElemSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  ElemDataRequests dataNeededByKernels_;
  stk::mesh::EntityRank entityRank_;
  unsigned nodesPerEntity_;
  int rhsSize_;
  const bool interleaveMEViews_;
};

} // namespace nalu
} // namespace Sierra

#endif

