/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <AveragingInfo.h>
#include <NaluParsing.h>

// basic c++
#include <stdexcept>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AveragingInfo - holder for averaging information held at TurbAvePP
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AveragingInfo::AveragingInfo() 
: computeReynoldsStress_(false),
  computeTke_(false),
  computeFavreStress_(false),
  computeFavreTke_(false),
  computeVorticity_(false),
  computeQcriterion_(false),
  computeLambdaCI_(false),
  computeMeanResolvedKe_(false)
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AveragingInfo::~AveragingInfo()
{
  // nothing to do
}


} // namespace nalu
} // namespace Sierra
