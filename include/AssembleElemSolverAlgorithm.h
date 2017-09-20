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
    const stk::topology &theTopo,
    bool interleaveMeViews = true);
  virtual ~AssembleElemSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  // topo and master element for this instance
  stk::topology topo_;
  MasterElement *meSCS_;
  MasterElement *meSCV_;

  ElemDataRequests dataNeededByKernels_;
  int rhsSize_;
  const bool interleaveMEViews_;
};

} // namespace nalu
} // namespace Sierra

#endif

