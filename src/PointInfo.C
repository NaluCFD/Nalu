/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <PointInfo.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{


//==========================================================================
// Class Definition
//==========================================================================
// PointInfo - contains point -> donor elements mapping
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
PointInfo::PointInfo(
  boundingPoint bPoint,
  const uint64_t localPointId,
  Point &ipCoords,
  Point &pointCoords,
  const int nDim)
  : bPoint_(bPoint),
    localPointId_(localPointId),
    ipCoordinates_(ipCoords),
    pointCoordinates_(pointCoords),
    owningElement_(),
    bestX_(1.0e16),
    bestXRef_(1.0e16),
    elemIsGhosted_(0),
    meSCS_(NULL)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
PointInfo::~PointInfo()
{
  // nothing to delete
}

} // namespace nalu
} // namespace sierra
