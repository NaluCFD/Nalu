/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
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
  const Realm &realm) 
  : AlgorithmDriver(realm)
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
SolverAlgorithmDriver::~SolverAlgorithmDriver()
{
  std::map<AlgorithmType, SolverAlgorithm *>::iterator ii;
  for( ii=solverAlgMap_.begin(); ii!=solverAlgMap_.end(); ++ii ) {
    Algorithm *theAlg = ii->second;
    delete theAlg;
  }
  
  std::map<AlgorithmType, SolverAlgorithm *>::iterator iid;
  for( iid=solverDirichAlgMap_.begin(); iid!=solverDirichAlgMap_.end(); ++iid ) {
    Algorithm *theAlg = iid->second;
    delete theAlg;
  }
  
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
SolverAlgorithmDriver::initialize_connectivity()
{

  std::map<AlgorithmType, SolverAlgorithm *>::iterator it;
  for ( it = solverAlgMap_.begin(); it != solverAlgMap_.end(); ++it ) {
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
  
  // assemble all interior and boundary contributions
  std::map<AlgorithmType, SolverAlgorithm *>::iterator it;
  for ( it = solverAlgMap_.begin(); it != solverAlgMap_.end(); ++it ) {
    it->second->execute();
  }
  
  // handle dirichlet
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd;
  for ( itd = solverDirichAlgMap_.begin(); itd != solverDirichAlgMap_.end(); ++itd ) {
    itd->second->execute();
  }

  post_work();
  
}


} // namespace nalu
} // namespace Sierra
