/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <AlgorithmDriver.h>

#include <Algorithm.h>
#include <Enums.h>

namespace sierra{
namespace nalu{

class Realm;

//==========================================================================
// Class Definition
//==========================================================================
// AlgorithmDriver - Drives algorithms
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AlgorithmDriver::AlgorithmDriver(
  Realm &realm)
  : realm_(realm)
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AlgorithmDriver::~AlgorithmDriver()
{
  std::map<AlgorithmType, Algorithm * >::iterator ii;
  for( ii=algMap_.begin(); ii!=algMap_.end(); ++ii ) {
    Algorithm *theAlg = ii->second;
    delete theAlg;
  }
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AlgorithmDriver::execute()
{
  pre_work();

  // assemble
  std::map<AlgorithmType, Algorithm *>::iterator it;
  for ( it = algMap_.begin(); it != algMap_.end(); ++it ) {
    it->second->execute();
  }

  post_work();

}


} // namespace nalu
} // namespace Sierra
