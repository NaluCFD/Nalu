/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <SolverAlgorithmDriver.h>

#include <AlgorithmDriver.h>
#include <Enums.h>
#include <SolverAlgorithm.h>

namespace sierra{
namespace nalu{

class Realm;

//==========================================================================
// Class Definition
//==========================================================================
// SolverAlgorithmDriver - Drives algorithms
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SolverAlgorithmDriver::SolverAlgorithmDriver(
  Realm &realm) 
  : AlgorithmDriver(realm)
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
SolverAlgorithmDriver::~SolverAlgorithmDriver()
{
  std::map<std::string, SolverAlgorithm *>::iterator is;
  for( is=solverAlgorithmMap_.begin(); is != solverAlgorithmMap_.end(); ++is ) {
    Algorithm *theAlg = is->second;
    delete theAlg;
  }

  std::map<AlgorithmType, SolverAlgorithm *>::iterator ii;
  for( ii=solverAlgMap_.begin(); ii!=solverAlgMap_.end(); ++ii ) {
    Algorithm *theAlg = ii->second;
    delete theAlg;
  }

  for( ii=solverConstraintAlgMap_.begin(); ii!=solverConstraintAlgMap_.end(); ++ii ) {
    Algorithm *theAlg = ii->second;
    delete theAlg;
  }
  
  for( ii=solverDirichAlgMap_.begin(); ii!=solverDirichAlgMap_.end(); ++ii ) {
    Algorithm *theAlg = ii->second;
    delete theAlg;
  }
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
SolverAlgorithmDriver::initialize_connectivity()
{
  std::map<std::string, SolverAlgorithm *>::iterator itc;
  for ( itc = solverAlgorithmMap_.begin(); itc != solverAlgorithmMap_.end(); ++itc ) {
    itc->second->initialize_connectivity();
  }
  std::map<AlgorithmType, SolverAlgorithm *>::iterator it;
  for ( it = solverAlgMap_.begin(); it != solverAlgMap_.end(); ++it ) {
    it->second->initialize_connectivity();
  }
  for ( it = solverConstraintAlgMap_.begin(); it != solverConstraintAlgMap_.end(); ++it ) {
    it->second->initialize_connectivity();
  }

  for ( it = solverDirichAlgMap_.begin(); it != solverDirichAlgMap_.end(); ++it ) {
    it->second->initialize_connectivity();
  }

}

//--------------------------------------------------------------------------
//-------- pre_work---------------------------------------------------------
//--------------------------------------------------------------------------
void
SolverAlgorithmDriver::pre_work()
{
  // might set initial guess
}


//--------------------------------------------------------------------------
//-------- post_work--------------------------------------------------------
//--------------------------------------------------------------------------
void
SolverAlgorithmDriver::post_work()
{
  // might provide residual
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
SolverAlgorithmDriver::execute()
{
  pre_work();
  
  // assemble all interior and boundary contributions; consolidated homogeneous approach
  std::map<std::string, SolverAlgorithm *>::iterator itc;
  for ( itc = solverAlgorithmMap_.begin(); itc != solverAlgorithmMap_.end(); ++itc ) {
    itc->second->execute();
  }

  // assemble all interior and boundary contributions
  std::map<AlgorithmType, SolverAlgorithm *>::iterator it;
  for ( it = solverAlgMap_.begin(); it != solverAlgMap_.end(); ++it ) {
    it->second->execute();
  }
  
  // handle constraint (will zero out entire row and process constraint)
  for ( it = solverConstraintAlgMap_.begin(); it != solverConstraintAlgMap_.end(); ++it ) {
    it->second->execute();
  }

  // handle dirichlet
  for ( it = solverDirichAlgMap_.begin(); it != solverDirichAlgMap_.end(); ++it ) {
    it->second->execute();
  }

  post_work();
  
}


} // namespace nalu
} // namespace Sierra
