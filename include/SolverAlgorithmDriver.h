/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SolverAlgorithmDriver_h
#define SolverAlgorithmDriver_h

#include<AlgorithmDriver.h>
#include<Enums.h>

#include<map>

namespace sierra{
namespace nalu{

class Realm;
class SolverAlgorithm;

class SolverAlgorithmDriver : public AlgorithmDriver
{
public:

  SolverAlgorithmDriver(
    Realm &realm);
  virtual ~SolverAlgorithmDriver();

  virtual void initialize_connectivity();
  virtual void pre_work();
  virtual void execute();
  virtual void post_work();
  
  // different types of algorithms... interior/flux; constraints and dirichlet
  std::map<std::string, SolverAlgorithm *> solverAlgorithmMap_;
  std::map<AlgorithmType, SolverAlgorithm *> solverAlgMap_;
  std::map<AlgorithmType, SolverAlgorithm *> solverConstraintAlgMap_;
  std::map<AlgorithmType, SolverAlgorithm *> solverDirichAlgMap_;
};

} // namespace nalu
} // namespace Sierra

#endif
